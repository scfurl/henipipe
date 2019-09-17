#!/usr/bin/python

import sys
try:
    from collections import OrderedDict
except ImportError: #python 2.6 or 3.6+
    if sys.version_info >= (3,6):
        OrderedDict = dict
    else:
        from ordereddict import OrderedDict

import os
from itertools import groupby
from subprocess import Popen, PIPE
from io import TextIOWrapper
import re
import argparse
from six import PY3, string_types

try:
    from multiprocessing.dummy.connection import Connection
except ImportError: #python2
    from _multiprocessing import Connection

__version__ = '0.1'

class DefaultOrderedDict(OrderedDict):
    def __init__(self, default, items=[]):
        super(DefaultOrderedDict, self).__init__(items)
        self._default = default

    def __missing__(self, key):
        self[key] = value = self._default()
        return value


class GenomicOrder(object):
    def __gt__(self, other):
        if self.rname != other.rname:
            return self.rname > other.rname
        return self.pos > other.pos

    def __lt__(self, other):
        if self.rname != other.rname:
            return self.rname < other.rname
        return self.pos < other.pos

    def __eq__(self, other):
        return self.rname == other.rname and self.pos == other.pos


class Reader(object):
    """ Read SAM/BAM format file as an iterable. """
    def __init__(self, f, regions=False, kind=None, samtools_path="samtools"):
        ext = None
        self.samtools_path = samtools_path
        self.spool = None  # use this to catch alignment during reader scraping
        self.type = 'sam'
        try:
            self._f_name = f.name
            _, ext = os.path.splitext(f.name)
            if f.name == '<stdin>':  # stdin stream
                self._sam_init(f)
            elif (ext is not None and ext.lower()) == '.bam' or (kind is not None and kind.lower() == 'bam'):
                self._bam_init(f, regions)
                self.type = 'bam'
            elif (ext is not None and ext.lower()) == '.sam' or (kind is not None and kind.lower() == 'sam'):
                self._sam_init(f)
            else:
                self._sam_init(f)
            if (regions and (ext is not None and ext.lower() != '.bam') and kind is None) or (regions and kind is not None and kind.lower() != 'bam'):
                self.__exit__()
                raise ValueError("Region support requires bam file.")
        except AttributeError:
            self._f_name = None
            if isinstance(f, Connection):
                self._pipe_init(f)
            else:
                self._sam_init(f)

    def _pipe_init(self, f):
        header = []
        for line in iter(f.recv, ''):
            if line[0] == '@':
                header.append(line.rstrip('\n\r'))
            else:
                self.spool = line
                break
        self.header_as_dict(header)
        self.f = iter(f.recv, '')
        self._conn = 'pipe'

    def _sam_init(self, f):
        header = []
        self.f = f
        for line in self.f:
            if line[0] == '@':
                header.append(line.rstrip('\n\r'))
            else:
                self.spool = line
                break
        self.header_as_dict(header)
        self._conn = 'file'

    def _bam_init(self, f, regions):
        pline = [self.samtools_path, 'view', '-H', f.name]
        try:
            p = Popen(pline, bufsize=-1, stdout=PIPE,
                      stderr=PIPE)
        except OSError:
            raise OSError('Samtools must be installed for BAM file support!\n')
        self.header_as_dict([line.decode('utf-8').rstrip('\n\r') for line in p.stdout])
        p.wait()
        if regions:
            try:
                open(''.join([f.name, '.bai']))
            except EnvironmentError:
                sys.stderr.write("BAM index not found. Attempting to index file.\n")
                index_p = Popen([self.samtools_path, 'index', f.name], stdout=PIPE, stderr=PIPE)

                _, err = index_p.communicate()
                if index_p.returncode > 0 or re.search("fail", str(err)):
                    raise OSError("Indexing failed. Is the BAM file sorted?\n")
                else:
                    sys.stderr.write("Index created successfully.\n")
            pline = [self.samtools_path, 'view', f.name, regions]
        else:
            pline = [self.samtools_path, 'view', f.name]
        self.p = Popen(pline, bufsize=-1, stdout=PIPE,
                  stderr=PIPE)
        if PY3:
            self.f = TextIOWrapper(self.p.stdout)
        else:
            self.f = self.p.stdout

        self._conn = 'proc'

    def next(self):
        """ Returns the next :class:`.Sam` object """
        try:
            if self.spool:  # this will be the first alignment in a SAM file or stream
                line = self.spool.rstrip('\n\r')
                self.spool = None
            else:
                line = next(self.f).rstrip('\n\r')
            if line == '':
                raise StopIteration
            fields = line.split('\t')
            required = fields[:11]
            tags = fields[11:]
            return Sam(*required, tags=tags)
        except StopIteration:
            raise StopIteration

    def __next__(self):
        return self.next()

    def __iter__(self):
        return self

    def __len__(self):
        """ Returns the number of reads in an indexed BAM file.
        Not implemented for SAM files. """
        if self.type != 'bam':
            raise NotImplementedError("len(Reader) is only implemented for BAM files.")
        elif self.type == 'bam':
            return sum(bam_read_count(self._f_name, self.samtools_path))

    def subsample(self, n):
        """ Returns an interator that draws every nth read from
        the input file. Returns :class:`.Sam`. """
        for i, line in enumerate(self.f):
            if i % n == 0:
                fields = line.split('\t')
                required = fields[:11]
                tags = fields[11:]
                yield Sam(*required, tags=tags)

    def header_as_dict(self, header):
        """ Parse the header list and return a nested dictionary. """
        self.header = DefaultOrderedDict(OrderedDict)
        for line in header:
            line = line.split('\t')
            key, fields = (line[0], line[1:])
            try:
                self.header[key][fields[0]] = fields[1:]
            except IndexError:
                self.header[key][fields[0]] = ['']

    @property
    def seqs(self):
        """ Return just the sequence names from the @SQ library as a generator. """
        for key in self.header['@SQ'].keys():
            yield key.split(':')[1]

    def tile_genome(self, width):
        """ Return a generator of UCSC-style regions tiling ``width``. """
        assert isinstance(width, int)
        for k, v in self.header['@SQ'].items():
            rname = k.split(':')[1]
            seqlength = v[0].split(':')[1]
            for region in tile_region(rname, 1, int(seqlength), width):
                yield region

    def close(self):
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        if self._conn == 'file':
            self.f.close()
        if self._conn == 'proc':
            self.f.close()
            self.p.terminate()


class Writer(object):
    """ Write SAM/BAM format file from :class:`.Sam` objects. """
    def __init__(self, f, header=None):
        try:
            _, ext = os.path.splitext(f.name)
            if ext == '.bam':
                raise NotImplementedError('Bam writing support is not implemented.\n') # why not just pipe to samtools?
        except AttributeError:
            pass
        self.file = f
        if header is not None:
            self.header = DefaultOrderedDict(OrderedDict)
            self._merge_header(header)
        else:
            self.header = DefaultOrderedDict(OrderedDict)
            self.header['@HD']['VN:1.0'] = ['SO:unknown']
        self._header_dict_format()

    def _merge_header(self, header):
        for key, values in header.items():
            for k, v in values.items():
                self.header[key][k] = v

    def _header_dict_format(self):
        for key, value in self.header.items():
            for k, v in value.items():
                tags = '\t'.join(v)
                self.file.write('{key}\t{k}\t{tags}\n'.format(**locals()))

    def write(self, sam):
        """ Write the string representation of the ``sam`` :class:`.Sam` object. """
        self.file.write(str(sam))

    def close(self):
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


class Sam(GenomicOrder):
    """ Object representation of a SAM entry. """
    # https://github.com/samtools/hts-specs/blob/da805be01e2ceaaa69fdde9f33c5377bf9ee6369/SAMv1.tex#L383
    # operations that consume the reference
    _cigar_ref = set(('M', 'D', 'N', '=', 'X', 'EQ'))
    # operations that consume the query
    _cigar_query = set(('M', 'I', 'S', '=', 'X', 'EQ'))
    # operations that do not represent an alignment
    _cigar_no_align = set(('H', 'P'))
    _valid_cigar = _cigar_ref | _cigar_query | _cigar_no_align
    # operations that can be represented as aligned to the reference
    _cigar_align = _cigar_ref & _cigar_query
    # operations that only consume the reference
    _cigar_ref_only = _cigar_ref - _cigar_align
    # operations that only consume the query
    _cigar_query_only = _cigar_query - _cigar_align

    def __init__(self, qname='', flag=4, rname='*', pos=0, mapq=255, cigar='*', rnext='*', pnext=0, tlen=0, seq='*', qual='*', tags=[]):
        self.qname = qname
        self.flag = int(flag)
        self.rname = rname
        self.pos = int(pos)
        self.mapq = int(mapq)
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = int(pnext)
        self.tlen = int(tlen)
        self.seq = seq
        self.qual = qual
        self._tags = tags
        self._cache = dict()

    def __str__(self):
        """ Returns the string representation of a SAM entry. Correspondes to one line
        in the on-disk format of a SAM file. """
        if self.tags:
            tag_fields = '\t'.join([encode_tag(tag, self.tags[tag]) for tag in sorted(self.tags.keys())])
        else:
            tag_fields = '\t'.join(self._tags)
        return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n'.format(self.qname,
                                                                                     str(self.flag),
                                                                                     self.rname,
                                                                                     str(self.pos),
                                                                                     str(self.mapq),
                                                                                     self.cigar,
                                                                                     self.rnext,
                                                                                     str(self.pnext),
                                                                                     str(self.tlen),
                                                                                     self.seq,
                                                                                     self.qual,
                                                                                     tag_fields)
    def __repr__(self):
        return "Sam({0}:{1}:{2})".format(self.rname, self.pos, self.qname)

    def __len__(self):
        """ Returns the length of the portion of ``self.seq`` aligned to the reference. Unaligned reads will
        have len() == 0. Insertions (I) and soft-clipped portions (S) will not contribute to the aligned length. 

        >>> x = Sam(cigar='8M2I4M1D3M4S')
        >>> len(x)
        16
        """
        return sum(c[0] for c in self.cigars if c[1] in self._cigar_ref)

    def __getitem__(self, tag):
        """ Retreives the SAM tag named "tag" as a tuple: (tag_name, data). The
        data type of the tag is interpreted as the proper Python object type.

        >>> x = Sam(tags=['NM:i:0', 'ZZ:Z:xyz'])
        >>> x['NM']
        0
        >>> x['ZZ']
        'xyz'
        """
        return self.tags[tag]

    def __setitem__(self, tag, data):
        """ Stores the SAM tag named "tag" with the value "data". The
        data type of the tag is interpreted from the Python object type.

        >>> x = Sam(tags=[])
        >>> x['NM'] = 0
        >>> x['NM']
        0
        """
        self.tags[tag] = data

    def index_of(self, pos):
        """ Return the relative index within the alignment from a genomic position 'pos' """
        i = pos - self.pos
        if i >= 0:
            return i
        else:
            raise IndexError("Position {0:n} not in {1}.".format(pos, self.qname))

    def get(self, key, default_value):
        try:
            return self[key]
        except KeyError:
            return default_value

    def cigar_split(self):
        # https://github.com/brentp/bwa-meth
        if self.cigar == "*":
            yield (0, None)
            raise StopIteration
        cig_iter = groupby(self.cigar, lambda c: c.isdigit())
        for _, n in cig_iter:
            op = int("".join(n)), "".join(next(cig_iter)[1])
            if op[1] in self._valid_cigar:
                yield op
            else:
                raise ValueError("CIGAR operation %s in record %s is invalid." % (op[1], self.qname))

    def gapped(self, attr, gap_char='-'):
        """ Return a :class:`.Sam` sequence attribute or tag with all
        deletions in the reference sequence represented as 'gap_char' and all
        insertions in the reference sequence removed. A sequence could
        be :class:``Sam.seq``, ``Sam.qual``, or any :class:`.Sam` tag that
        represents an aligned sequence, such as a methylation tag for bisulfite
        sequencing libraries.

        >>> x = Sam(*'r001\t99\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\t*'.split())
        >>> x.gapped('seq')
        'TTAGATAAGATA-CTG'
        >>> x = Sam(*'r001\t99\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\t*'.split(), tags=['ZM:Z:.........M....M.M'])
        >>> x.gapped('ZM')
        '............-M.M'
        """
        try:
            ungapped = getattr(self, attr)
        except AttributeError:
            ungapped = self[attr]  # get dictionary key (tag) if attribute is missing
        if len(ungapped) != len(self.seq):
            raise ValueError("The length of the '%s' attribute is not equal to the length of Sam.seq!" % attr) 
        gapped = []
        i = 0
        for n, t in self.cigars:
            if t in self._cigar_align:
                gapped.extend(ungapped[i:i + n])
                i += n
            elif t in self._cigar_ref_only:
                gapped.extend([gap_char] * n)
            elif t in self._cigar_query_only:
                i += n
            elif t in self._cigar_no_align:
                pass
        return ''.join(gapped)

    def parse_md(self):
        """ Return the ungapped reference sequence from the MD tag, if present.
        """
        try:
            return self._cache['parse_md']
        except KeyError:
            pass
        try:
            md = self['MD']
        except KeyError:
            raise KeyError('MD tag not found in SAM record.')
        ref_seq = list(self.gapped('seq'))
        md_match = re.findall(r"([0-9]+)\^?([A-Z]+)?", md)
        ref_seq_i = 0
        for i, b in md_match:
            ref_seq_i += int(i)
            for mismatch in b:
                try:
                    ref_seq[ref_seq_i] = mismatch
                except IndexError:
                    raise IndexError(locals())
                ref_seq_i += 1
        self._cache['parse_md'] = ref_seq
        return ref_seq

    @property
    def cigars(self):
        """ Returns the CIGAR string as a tuple.

        >>> x = Sam(cigar='8M2I4M1D3M')
        >>> x.cigars
        ((8, 'M'), (2, 'I'), (4, 'M'), (1, 'D'), (3, 'M'))
        """
        try:
            return self._cache['cigars']
        except KeyError:
            self._cache['cigars'] = tuple(self.cigar_split())
            return self._cache['cigars']

    @property
    def tags(self):
        """ Parses the tags string to a dictionary if necessary.

        >>> x = Sam(tags=['XU:Z:cgttttaa', 'XB:Z:cttacgttaagagttaac', 'MD:Z:75', 'NM:i:0', 'NH:i:1', 'RG:Z:1'])
        >>> sorted(x.tags.items(), key=lambda x: x[0])
        [('MD', '75'), ('NH', 1), ('NM', 0), ('RG', '1'), ('XB', 'cttacgttaagagttaac'), ('XU', 'cgttttaa')]
        """
        try:
            return self._cache['tags']
        except KeyError:
            self._cache['tags'] = parse_sam_tags(self._tags)
            return self._cache['tags']

    @property
    def paired(self):
        """ Returns True if the read is paired and
        each segment properly aligned according to the aligner. """
        return bool(self.flag & 0x2)

    @property
    def mapped(self):
        """ Returns True of the read is mapped. """
        return not (self.flag & 0x4)

    @property
    def secondary(self):
        """ Returns True if the read alignment is secondary. """
        return bool(self.flag & 0x100)

    @property
    def reverse(self):
        """ Returns True if ``Sam.seq`` is being reverse complemented. """
        return bool(self.flag & 0x10)

    @property
    def passing(self):
        """ Returns True if the read is passing filters, such as platform/vendor quality controls. """
        return not bool(self.flag & 0x200)

    @property
    def duplicate(self):
        """ Returns True if the read is a PCR or optical duplicate. """
        return bool(self.flag & 0x400)

    @property
    def coords(self):
        """ Return a list of genomic coordinates for the gapped alignment. """
        return range(self.pos, self.pos + len(self))

    @property
    def safename(self):
        """Return ``Sam.qname`` without paired-end identifier if it exists"""
        if self.qname[-2] == '/':
            return self.qname[:-2]
        else:
            return self.qname


def parse_sam_tags(tagfields):
    """ Return a dictionary containing the tags """
    return dict([(tag, data) for tag, dtype, data in [decode_tag(x) for x in tagfields]])


def encode_tag(tag, data):
    """ Write a SAM tag in the format ``TAG:TYPE:data``. Infers the data type
    from the Python object type.

    >>> encode_tag('YM', '#""9O"1@!J')
    'YM:Z:#""9O"1@!J'
    """
    if isinstance(data, string_types):
        data_type = 'Z'
    elif isinstance(data, int):
        data_type = 'i'
    elif isinstance(data, float):
        data_type = 'f'
    else:
        raise NotImplementedError("Data {0} cannot be encoded as string, integer, or float tag.".format(data))
    value = ':'.join((tag, data_type, str(data)))
    return value


def decode_tag(tag_string):
    """ Parse a SAM format tag to a (tag, type, data) tuple. Python object
    types for data are set using the type code. Supported type codes are: A, i, f, Z, H, B

    >>> decode_tag('YM:Z:#""9O"1@!J')
    ('YM', 'Z', '#""9O"1@!J')
    >>> decode_tag('XS:i:5')
    ('XS', 'i', 5)
    >>> decode_tag('XF:f:100.5')
    ('XF', 'f', 100.5)
    """
    try:
        tag, data_type, data = tag_string.split(':')
    except ValueError:
        match = re.match(r'([A-Z]{2}):([iZfHB]):(\S+)', tag_string)
        tag = match.group(1)
        data_type = match.group(2)
        data = match.group(3)
    if data_type == 'i':
        return (tag, data_type, int(data))
    elif data_type == 'Z':
        return (tag, data_type, data)
    elif data_type == 'f':
        return (tag, data_type, float(data))
    elif data_type == 'A':  # this is just a special case of a character
        return (tag, data_type, data)
    elif data_type == 'H':
        raise NotImplementedError("Hex array SAM tags are currently not parsed.")
    elif data_type == 'B':
        raise NotImplementedError("Byte array SAM tags are currently not parsed.")
    else:
        raise NotImplementedError("Tag {0} cannot be parsed.".format(tag_string))


def tile_region(rname, start, end, step):
    """ Make non-overlapping tiled windows from the specified region in
    the UCSC-style string format.

    >>> list(tile_region('chr1', 1, 250, 100))
    ['chr1:1-100', 'chr1:101-200', 'chr1:201-250']
    >>> list(tile_region('chr1', 1, 200, 100))
    ['chr1:1-100', 'chr1:101-200']
    """
    while start + step <= end:
        yield '%s:%d-%d' % (rname, start, start + step - 1)
        start += step
    if start < end:
        yield '%s:%d-%d' % (rname, start, end)

def bam_read_count(bamfile, samtools_path="samtools"):
    """ Return a tuple of the number of mapped and unmapped reads in a BAM file """
    p = Popen([samtools_path, 'idxstats', bamfile], stdout=PIPE)
    mapped = 0
    unmapped = 0
    for line in p.stdout:
        rname, rlen, nm, nu = line.rstrip().split()
        mapped += int(nm)
        unmapped += int(nu)
    return (mapped, unmapped)


class fragment(object):
    def __init__(self, loc, bed_start, end, fraglen, R1_mapq, R2_mapq, R1_flag, R2_flag, direction):
        self.loc = loc
        self.bed_start = bed_start
        self.end = end
        self.fraglen = fraglen
        self.R1_mapq = R1_mapq
        self.R2_mapq = R2_mapq
        self.R1_flag = R1_flag
        self.R2_flag = R2_flag
        self.direction = direction
    def __repr__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (self.loc,self.bed_start,self.end,self.fraglen,self.R1_mapq,self.R2_mapq,self.R1_flag, self.R2_flag, self.direction)
    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (self.loc,self.bed_start,self.end,self.fraglen,self.R1_mapq,self.R2_mapq,self.R1_flag, self.R2_flag, self.direction)

class fragments(object):
    def __init__(self, in_file, out_file=None, skip_dups = True, filter_threshes = None):
        #self.in_file = in_file
        self.read_count=0
        self.fragment_count = 0
        if filter_threshes is None:
            self.filter_threshes = [float("-inf"),float("inf")]
        else:
            self.filter_threshes = filter_threshes
        self.passed_filter = 0
        if skip_dups:
            self.forward_list = [99, 163, 355, 419]
            self.plus_list = [99, 355]
        else:
            self.forward_list = [99, 163, 355, 419, 1123, 1187]
            self.plus_list = [99, 355, 1123]
        if out_file is None:
            self.return_fragments(in_file)
        else:
            self.write_fragments(in_file, out_file)
        # self.generator = self.parse_sam_file(self.in_sam)

    def load_fragments(self, in_file):

        in_sam = Reader(in_file)
        iterator = iter(in_sam)
        done_looping = False
        fragments = []
        while not done_looping:
            try:
                read1 = next(iterator)
                self.read_count += 1
                read2 = next(iterator)
                self.read_count += 1
            except StopIteration:
                done_looping = True
                #print("Processed %s reads.  Found %s usable fragments" % (self.read_count, self.fragment_count))
            else:
                return_val = self.qualify_reads(read1, read2)
                if return_val is not None:
                    fragments.append(return_val)
                #print(str(self.qualify_reads(read1, read2)))
        return fragments

    def return_fragments(self, in_file):
        in_sam = Reader(in_file)
        iterator = iter(in_sam)
        done_looping = False
        while not done_looping:
            try:
                read1 = next(iterator)
                self.read_count += 1
                read2 = next(iterator)
                self.read_count += 1
            except StopIteration:
                done_looping = True
            except IOError:
                done_looping = True
                #print("Processed %s reads.  Found %s usable fragments" % (self.read_count, self.fragment_count))
            else:
                return_val = self.qualify_reads(read1, read2)
                if return_val is not None:
                    sys.stdout.write(str(return_val))
                #print(str(self.qualify_reads(read1, read2)))

    def write_fragments(self, in_file, out_file):

        in_sam = Reader(in_file)
        iterator = iter(in_sam)
        done_looping = False
        fragments = []
        while not done_looping:
            try:
                read1 = next(iterator)
                self.read_count += 1
                read2 = next(iterator)
                self.read_count += 1
            except StopIteration:
                done_looping = True
                print("\n[SAM2BED] Output: \n Processed %s read pairs.  Found %s usable fragments.\nOf these, %s passed_filter.\n" % (self.read_count/2, self.fragment_count, self.passed_filter))
            else:
                return_val = self.qualify_reads(read1, read2, filter = self.filter_threshes)
                if return_val is not None:
                    out_file.writelines(str(line) for line in [return_val])

    # def qualify_reads_old(self, read1, read2, filter = [float("-inf"),float("inf")]):
    #     direction = "Unk"
    #     r1_unique = r2_unique = False
    #     if read1.flag in self.forward_list:
    #         if "XS" not in read1.tags:
    #             r1_unique=True
    #             if "XS" not in read2.tags:
    #                 r2_unique=True
    #                 begin = read1.pos
    #                 end = (read2.pos + len(read2.seq))
    #                 direction = "+"
    #     if read2.flag in self.forward_list:
    #         if "XS" not in read2.tags:
    #             r2_unique=True
    #             if "XS" not in read1.tags:
    #                 r1_unique=True
    #                 begin = read2.pos
    #                 end = read1.pos - len(read1.seq)
    #                 direction = "-"
    #     if r1_unique and r2_unique and direction is not "Unk" and (end - begin) >0 :
    #         self.fragment_count+=1
    #         if (end - begin) > filter[0] and (end - begin) < filter[1]:
    #             self.passed_filter+=1
    #             return fragment(read1.rname, begin, end, end - begin, read1.mapq, read2.mapq, read1.flag, read2.flag, direction)

    # Sam fields: qname='', flag=4, rname='*', pos=0, mapq=255, cigar='*', rnext='*', pnext=0, tlen=0, seq='*', qual='*', tags=[]
    def qualify_reads(self, read1, read2, filter = [float("-inf"),float("inf")]):
        direction = "Unk"
        r1_unique = r2_unique = False
        if read1.flag in self.forward_list:
            if "XS" not in read1.tags:
                r1_unique=True
                if "XS" not in read2.tags:
                    r2_unique=True
                    begin = read1.pos
                    end = read1.pos + read1.tlen
                    direction = "+"
        if read2.flag in self.forward_list:
            if "XS" not in read2.tags:
                r2_unique=True
                if "XS" not in read1.tags:
                    r1_unique=True
                    begin = read2.pos
                    end = read2.pos + read2.tlen
                    direction = "-"
        if r1_unique and r2_unique and direction is not "Unk" and (end - begin) >0 and read1.rnext == "=" and read2.rnext == "=":
            self.fragment_count+=1
            if (end - begin) > filter[0] and (end - begin) < filter[1]:
                self.passed_filter+=1
                return fragment(read1.rname, (begin -1) , (end - 1) , end - begin, read1.mapq, read2.mapq, read1.flag, read2.flag, direction)


def run_sam2bed():
    parser = argparse.ArgumentParser('A script for converting sam file to bed for eventual bedgraph conversion')
    parser.add_argument('sam_bam', type=str, default = sys.stdin, help='sam or bam fil or stdin')
    parser.add_argument('--out_bed', '-o', default = None, help='bed output', required=False)
    parser.add_argument('--filter_high', '-fh', type = int, default = None, help='filter_high', required=False)
    parser.add_argument('--filter_low', '-fl', default = None, type = int,  help='filter_high', required=False)
    args = parser.parse_args()
    if args.sam_bam=="-":
        args.sam_bam = sys.stdin
    else:
        args.sam_bam = open(args.sam_bam, 'r')

    filter_threshes = []
    if args.filter_low is not None:
        filter_threshes.append(args.filter_low)
    else:
        filter_threshes.append(float("-inf"))
    if args.filter_high is not None:
        filter_threshes.append(args.filter_high)
    else:
        filter_threshes.append(float("inf"))


    if args.out_bed is not None:
        of = open(args.out_bed, "w")
        fragments(args.sam_bam, out_file=of, filter_threshes = filter_threshes)
        of.close()
    else:
        fragments(args.sam_bam, filter_threshes = filter_threshes)

if __name__ == '__main__':
    run_sam2bed()




