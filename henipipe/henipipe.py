#!/usr/bin/python
# PBS/SLURM cluster job submission in Python
#v 0.1

# A wrapper that performs the necessary pipeline for generation of CNR/CNT data
# Written by Scott Furlan with code inspiration from Andrew Hill's cellwrapper; 
# uses a custom python script sam2bed2.py which takes code from a fantastic 
# sam reader "simplesam" - https://github.com/mdshw5/simplesam

# runsheet columns
# 1 - 'sample' name of the sample REQUIRED
# 2 - 'fasta' location of the Bowtie2 indexed fasta file REQUIRED
# 3 - 'spikein_fasta' location of the Bowtie2 indexed fasta file for spike_in normalization OPTIONAL
# 4 - 'fastq1' a tab seperated string of filenames denoting location of all R1 files for a sample REQUIRED
# 5 - 'fastq2' a tab seperated string of filenames denoting location of all R2 files for a sample REQUIRED
# 6 - 'bed_out' name of the location for the aligned and sorted bam file REQUIRED
# 7 - 'spikein_bed_out' name of the location for the aligned and sorted bam file OPTIONAL
# 8 - 'genome_sizes' REQUIRED
# 9 - 'bedgraph' file name of normalized bedgraph REQUIRED
# 10 - 'SEACR_key' sample key corresponding to sample groups to be run against an IgG (or other) contol.  all samples to be run against a control are given the same name and the control is labeleled with the an additional flag +'_CONTROL' (i.e. 4JS_CONTROL) OPTIONAL
# 11 - 'SEACR_out' file name of SEACR output OPTIONAL

import os
from subprocess import Popen, PIPE
import argparse
import time
import sys
import csv
import re
import string
import random
from itertools import chain, compress
#import pandas as pd


class SampleFactory:
    def __init__(self, *args, **kwargs):
        self.user = kwargs.get('user')
        self.cluster = kwargs.get('cluster')
        self.runsheet_data = kwargs.get('runsheet_data')
        self.debug = kwargs.get('debug')
        self.log_name = kwargs.get('log')
    def __call__():
        pass
        # if self.debug == False:
        #   open(self.log, 'w')

    def id_generator(self, size=6, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

    def generate_job(self):
        job_string=[]
        torun = len(self.runsheet_data)
        for i in range(torun):
            log_file = self.id_generator()
            job_name = self.job + "_" + self.runsheet_data[i]['sample']
            command = self.command[i]
            if self.cluster=="PBS":
                #to_append = '#!/bin/bash\n#PBS -N %s\n#PBS -l %s\n#PBS -j oe\n#PBS -o $PBS_O_WORKDIR/%s\n#PBS -A %s\ncd $PBS_O_WORKDIR\n%s' % (job_name, self.processor_line, self.log_name, self.user, command)
                #to_append = '#!/bin/bash\n#PBS -N %s\n#PBS -l %s\n#PBS -j oe\n#PBS -o $PBS_O_WORKDIR/%s\n#PBS -A %s\ncd $PBS_O_WORKDIR\n%s' % (job_name, self.processor_line, log_file, self.user, command)
                #to_append = "#!/bin/bash\n#PBS -N %s\n#PBS -l  %s\n#PBS -j oe\n#PBS -o $PBS_JOBDIR/%s\n#PBS -A %s\ncd $PBS_O_WORKDIR\n%s\nsed -e 's/^/[HENIPIPE] %s: /' $PBS_JOBDIR/%s >> %s\n" % (job_name, self.processor_line, log_file, self.user, command, job_name, log_file, self.log_name)
                to_append = "#!/bin/bash\n#PBS -N %s\n#PBS -l %s\n#PBS -j oe\n#PBS -o $PBS_O_WORKDIR/logtmp\n#PBS -A %s\ncd $PBS_O_WORKDIR\n{%s} 2>&1 | tee %s\nsed -e 's/^/[HENIPIPE] JOB: %s:\t\t/' %s >> %s\nrm %s\n" % (job_name, self.processor_line, self.user, command, log_file, job_name, log_file, self.log_name, log_file)
            if self.cluster=="SLURM":
                to_append = '#!/bin/bash\n#SBATCH --job-name=%s\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=1\n#SBATCH --mem-per-cpu=8000\n%s' % (job_name, command)
            job_string.append(to_append)
        return job_string

    def run_job(self):
        for script in self.script:
            if self.cluster=="PBS":
                if self.debug==False:
                    # Open a pipe to the command.
                    proc = Popen('qsub', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
                    if (sys.version_info > (3, 0)):
                        proc.stdin.write(script.encode('utf-8'))
                        out, err = proc.communicate()
                    else:
                        proc.stdin.write(script)
                        out, err = proc.communicate()
                # Print your job and the system response to the screen as it's submitted
                print(script)
                if self.debug==False:
                    print(out)
                    time.sleep(0.1)
            if self.cluster=="SLURM":
                if self.debug==False:
                    #open("temp.sh", 'w').write(script+'\n')
                    proc = Popen(['sbatch'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
                    out, err = proc.communicate(script)
                print(script)
                if self.debug==False:
                    print(out)
                    time.sleep(0.1)


class Align(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(Align, self).__init__(*args, **kwargs)
        self.bowtie_flags = kwargs.get('bowtie_flags')
        self.job = "HENIPIPE_ALIGN"
        self.filter_string = self.make_filter_string(kwargs.get('filter')[0], kwargs.get('filter')[1])
        self.norm_method = kwargs.get('norm_method')
        self.processor_line = self.align_processor_line()
        self.command = self.align_executable()
        self.script = self.generate_job()
    def __call__():
        pass

    def make_filter_string(self, low, high):
        if low is None and high is None:
            return("")
        if low is not None and high is None:
            return("-fl %s" % low)
        if low is None and high is not None:
            return("-fh %s" % high)
        if low is not None and high is not None:
            return("-fl %s -fh %s" % (low, highsCYB))

    def align_executable(self):
        commandline=""
        command = []
        for sample in self.runsheet_data:
            fastq1=re.sub('\t', ',', sample['fastq1'])
            fastq2=re.sub('\t', ',', sample['fastq2'])
            JOBSTRING = self.id_generator(size=10)
            sam2bed_string = """| sam2bed - -o %s %s""" % (sample['bed_out']+'tmp', self.filter_string)
            if self.cluster=="SLURM":
                modules = """\nsource /app/Lmod/lmod/lmod/init/bash\nml SEACR\nml bowtie2\nmodule load samtools\n"""
            else:
                modules = """\nmodule load python\nmodule load bowtie2\nmodule load samtools\necho '\nRunning Bowtie piped to sam2bed...\n[BOWTIE] Output:\n'\n"""
            norm_bowtie_flags='--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant -q --phred33 -I 10 -X 700'
            commandline = """bowtie2 %s -p 4 -1 %s -2 %s -x %s %s\n""" % (self.bowtie_flags, fastq1, fastq2, sample['fasta'], sam2bed_string)
            commandline = commandline + """\necho 'Sorting Bed...\n'\nsort -k1,1 -k2n,2n %s > %s\n""" % (sample['bed_out']+'tmp', sample['bed_out'])
            commandline = commandline + """rm %s \n""" % (sample['bed_out']+'tmp')
            if self.norm_method == "spike_in":
                commandline = commandline + """echo '\n[BOWTIE] Running Bowtie piped to sam2bed.py for spikein... Output:\n'\nbowtie2 %s -p 4 -1 %s -2 %s -x %s | /home/sfurla/Scripts/sam2bed.py - -o %s\n""" % (norm_bowtie_flags, fastq1, fastq2, sample['spikein_fasta'], sample['spikein_bed_out']+'tmp')
                commandline = commandline + """\necho 'Sorting Bed for spikein...\n'sort -k1,1 -k2n,2n %s > %s\n""" % (sample['spikein_bed_out']+'tmp', sample['spikein_bed_out'])
                commandline = commandline + """rm %s \n""" % (sample['spikein_bed_out']+'tmp')
            commandline = modules + commandline
            command.append(commandline)
        return command


    def align_processor_line(self):
        if self.cluster=="PBS":
            return """select=1:mem=16GB:ncpus=4"""
        if self.cluster=="SLURM":
            return ''


class Norm(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(Norm, self).__init__(*args, **kwargs)
        self.job = "HENIPIPE_NORM"
        norm_method = kwargs.get('norm_method')
        self.processor_line = self.norm_processor_line()
        self.norm_values = self.get_norm_values(method = norm_method)
        self.command = self.norm_executable()
        self.script = self.generate_job()
    def __call__():
        pass

    # def get_norm_values(self, method):
    #   if method=="read_count":
    #       for sample in self.runsheet_data:
    #           count = 0
    #           for line in open(sample['bed_out']).xreadlines(): count += 1
    #           ncount = 0
    #           for line in open(sample['spikein_bed_out']).xreadlines(): ncount += 1
    #           sample['scale_factor'] = (1/(float(count)/float(ncount))*100)
    #   if method=="coverage":
    #       for sample in self.runsheet_data:
    #           count = 0
    #           with open(sample['bed_out']) as infile:
    #               for line in infile:
    #                   count += int(line.split("\t")[4])
    #           ncount = 0
    #           with open(sample['spikein_bed_out']) as infile:
    #               for line in infile:
    #                   ncount += int(line.split("\t")[4])
    #           sample['scale_factor'] = (1/(float(count)/float(ncount))*100)# print(sample_norm_values)
    #   return

    def get_norm_values(self, method):
        if method=="read_count":
            for sample in self.runsheet_data:
                count = 0
                for line in open(sample['bed_out']).xreadlines(): count += 1
                sample['scale_factor'] = (1/float(count))*float(10**7)
        if method=="coverage":
            for sample in self.runsheet_data:
                count = 0
                with open(sample['bed_out']) as infile:
                    for line in infile:
                        count += int(line.split("\t")[4])
                sample['scale_factor'] = (1/float(count))*float(10**10)
        if method=="spike_in":
            for sample in self.runsheet_data:
                count = 0
                for line in open(sample['bed_out']).xreadlines(): count += 1
                ncount = 0
                for line in open(sample['spikein_bed_out']).xreadlines(): ncount += 1
                sample['scale_factor'] = ((float(count)/float(ncount))/100)
        return


    def norm_executable(self):
        commandline=""
        command = []
        for sample in self.runsheet_data:
            JOBSTRING = self.id_generator(size=10)
            if self.cluster=="SLURM":
                modules = """\nsource /app/Lmod/lmod/lmod/init/bash\n"""
            else:
                modules = """\nmodule load bedtools\n"""
                #modules = """\nmodule load bedtools\necho '\n[GENOMECOVERAGEBED] Making bedgraph... Output (None is good):\n'"""
            commandline = """genomeCoverageBed -bg -i %s -g %s -scale %s -trackline | pyWriter - %s\n""" % (sample['bed_out'], sample['genome_sizes'], sample['scale_factor'], sample['bedgraph'])
            commandline = modules + commandline
            command.append(commandline)
        return command


    def norm_processor_line(self):
        if self.cluster=="PBS":
            return """select=1:mem=8GB:ncpus=1"""
        if self.cluster=="SLURM":
            return ''

class SEACR(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(SEACR, self).__init__(*args, **kwargs)
        self.job = "HENIPIPE_SEACR"
        self.method = kwargs.get('stringency')
        self.norm = kwargs.get('norm')
        self.runsheet_data = self.SEACR_match(pare_down = kwargs.get('pare_down'))
        self.processor_line = self.SEACR_processor_line()
        self.command = self.SEACR_executable()
        self.script = self.generate_job()
    def __call__():
        pass

    def SEACR_match(self, pare_down):
        #will need to change this when multiple selections are implemented; for now just allow user to specify sample, then find control
        desired_samples = [self.runsheet_data[i].get("sample") for i in pare_down]
        sk = [i.get('SEACR_key') for i in self.runsheet_data]
        controls_b = [bool(re.search(r'._CONTROL$', i)) for i in sk]
        controls = list(compress(self.runsheet_data, controls_b))
        samples_b = [not i for i in controls_b]
        samples = list(compress(self.runsheet_data, samples_b))
        samples = [i for i in samples if i.get("sample") in desired_samples]
        for sample in samples:
            control_name = sample.get('SEACR_key')+"_CONTROL"
            control_bed = next(item for item in controls if item["SEACR_key"] == control_name).get('bedgraph')
            sample.update( {'SEACR_in' : sample.get('bedgraph')})
            sample.update( {'SEACR_control' : control_bed})
        return samples


    def SEACR_executable(self):
        commandline=""
        command = []
        #print("Runmode is " + self.runmode)
        for sample in self.runsheet_data:
            JOBSTRING = self.id_generator(size=10)
            if self.cluster=="SLURM":
                modules = """\nsource /app/Lmod/lmod/lmod/init/bash\nml SEACR\n"""
            else:
                modules = """\nmodule load SEACR\nmodule load R\nmodule load bedtools\n"""
            #commandline = """echo '\n[SEACR] Running SEACR... Output:\n'bash /home/sfurla/develop/SEACR/SEACR_1.1.sh %s %s %s %s %s""" % (sample['SEACR_in'], sample['SEACR_control'], self.norm, self.method, sample['SEACR_out'])
            commandline = """bash echo '\n[SEACR] Running SEACR... Output:\n'\nbash SEACR_1.1.sh %s %s %s %s %s\n""" % (sample['SEACR_in'], sample['SEACR_control'], self.norm, self.method, sample['SEACR_out'])
            commandline = modules + commandline
            command.append(commandline)
        return command


    def SEACR_processor_line(self):
        if self.cluster=="PBS":
            return """select=1:mem=8GB:ncpus=2"""
        if self.cluster=="SLURM":
            return ''


def convert_windows_newlines(file_name):
    """`
    Helper function for converting windows newlines in a file as a preprocessing step in case users make samplesheet
    in excel.

    Args:
    file_name (str): name of file to convert to unix newlines.

    """
    # Read in the file in entirety first
    with open(file_name) as input:
        file_text = input.read()

    # Replace newlines produced by excel
    file_text = file_text.replace('\r', '\n')

    # Write the new output back to the file
    with open(file_name, 'w') as output:
        output.write(file_text)


def reduce_concat(x, sep=""):
    return functools.reduce(lambda x, y: str(x) + sep + str(y), x)

# def paste(*lists, sep=" ", collapse=None):
#   result = map(lambda x: reduce_concat(x, sep=sep), zip(*lists))
#   if collapse is not None:
#       return reduce_concat(result, sep=collapse)
#   return list(result)
# def print_iterator(it):
#   for x in it:
#       print(x, end='\n')
#   print(' ')  # for new line

def find_fastq_mate(dir, sample_flag=None, full_name=True):
    fastqs=[]
    fastq1=[]
    fastq2=[]
    for file in os.listdir(dir):
        if file.endswith(".fastq.gz"):
            fastqs.extend([file])
            if "_R1_" in file:
                fastq1.extend([file])
            if "_R2_" in file:
                fastq2.extend([file])
    fastq1_mate=[]
    for fastq in fastq1:
        #check if present
        put_R2=re.sub('_R1_', '_R2_', fastq)
        try:
            fastq1_mate.extend([fastq2[fastq2.index(put_R2)]])
        except ValueError:
            raise ValueError("Could not find matching file: "+put_R2+" for: "+fastq+" in "+dir)
    if full_name:
            for i in range(len(fastq1)):
                fastq1[i] = os.path.join(dir, fastq1[i])
                fastq1_mate[i] = os.path.join(dir, fastq1_mate[i])
    keys=['directory_long', 'directory_short','fastq1', 'fastq2', 'has_fastq']
    #values=[dir, os.path.basename(dir), "\t".join(os.path.join(dir, fastq1)), "\t".join(os.path.join(dir, fastq1_mate)), len(fastq1)>0]
    values=[dir, os.path.basename(dir), "\t".join(fastq1), "\t".join(fastq1_mate), len(fastq1)>0]
    return(dict(zip(keys, values)))



def find_colnames(runsheet, header=True):
    """
    Helper function for getting headers for a runsheet.
    Args:
        runsheet (file): file object to runsheet.  CSV file with header is supported.
    Yields:
        Colnames of CSV
    """
    if header==True:
        with open(runsheet, 'r') as f:
            reader = csv.reader(f)
            return(next(reader))            # read header

def make_runsheet(folder, sample_flag = "samples", output=None, fasta=None, spikein_fasta=None, genome_sizes=None, furlan=True):
    #folder = '/active/furlan_s/Data/CNR/190801_CNRNotch/fastq/mini/fastq'
    if furlan:
        genome_sizes = '/active/furlan_s/Data/CNR/190801_CNRNotch/fastq/sizes.genome'
        fasta = '/active/furlan_s/refs/hg38/bowtie_hg38_p12'
        spikein_fasta = '/active/furlan_s/refs/ecoli/GCF_000005845.2_ASM584v2'
    if output is None:
        output = os.path.join(os.getcwd(), "HeniPipeOut")
    ddir=[x[0] for x in os.walk(folder)]
    dat=list(map(find_fastq_mate, ddir))
    good_dat = [i for i in dat if i.get('has_fastq') is True]
    good_dat = [i for i in good_dat if re.compile(r'.*'+sample_flag).search(i.get('directory_short'))]

    for i in good_dat:
        i.update({'sample': i.get('directory_short'), \
            'bed_out': os.path.join(output, i.get('directory_short')+".bed"), \
            'spikein_bed_out': os.path.join(output, i.get('directory_short')+"_spikein.bed"), \
            'bedgraph': os.path.join(output, i.get('directory_short')+".bedgraph"), \
            'SEACR_key': i.get('directory_short'), \
            'SEACR_out': os.path.join(output, i.get('directory_short')+"_SEACR.bedgraph")})
        if fasta is None:
            i.update({'fasta': 'ADD_BOWTIE_INDEX_HERE - i.e. path to the Bowtie2 indexed FASTA genome file'})
        else: i.update({'fasta': fasta})
        if spikein_fasta is None:
            i.update({'spikein_fasta': 'ADD_BOWTIE_SPIKEIN_INDEX_HERE - i.e. path to the Bowtie2 indexed spike-in FASTA genome file'})
        else: i.update({'spikein_fasta': spikein_fasta})
        if genome_sizes is None:
            i.update({'genome_sizes': 'ADD GENOME SIZES FILE HERE'})
        else: i.update({'genome_sizes':  genome_sizes})
    #print(good_dat)
    keys = good_dat[0].keys()
    LOGGER.info("Writing runsheet to - "+os.path.join(output, 'runsheet.csv')+" ...")
    with open(os.path.join(output, 'runsheet.csv'), 'wb') as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(good_dat)

def parse_runsheet(runsheet, header=True, colnames=None):
    """
    Helper function for parsing runsheets provided by user.
    Args:
        runsheet (file): file object to runsheet.
    Yields:
        Each line of the run as a dict by column name.
        Columns: as listed in csv
    """
    if header:
        columns = find_colnames(runsheet, header=header)
    if header==False and colnames is not None:
        columns = colnames
    if header==False and colnames is None:
        raise ValueError('No header indicated in runsheet and no colnames provided')

    entries = []
    with open(runsheet, 'r') as f:
        reader = csv.reader(f)
        if header==True:
            next(reader)            # skip header
        for line in f:
                #print(line)
            entries = [entry.strip() for entry in line.strip().split(',')]
            if len(entries) != len(columns):
                raise ValueError('Sample sheet does not match expected columns. Expects: %s' % ','.join(columns))
            entry_dict = dict(zip(columns, entries))
            yield entry_dict

def check_runsheet_parameter(runsheet, parameter, verbose=False):
    for i in runsheet:
        if i.get(parameter) is None:
            if verbose: print("no data for paramter "+parameter)
            return 0
        if i.get(parameter) is "":
            if vebose: print("header present, but no or incomplete data for "+parameter)
            return 0
    if verbose: print("runsheet_okay_for_parameter_"+parameter)
    return 1

def check_runsheet(args, runsheet, verbose=False):
    required_args = ['sample', 'fastq1', 'fastq2', 'bed_out', 'bedgraph', 'genome_sizes']
    if args.norm_method == "spike_in":
        required_args.extend(['spikein_fasta', 'spikein_bed_out'])
    if args.job == "SEACR":
        required_args.extend(['SEACR_key', 'SEACR_out'])
    has_data = []
    for i in required_args:
        has_data.append(check_runsheet_parameter(runsheet, i, verbose = verbose))
    has_data = dict(zip(required_args, has_data))
    listOfKeys = getKeysByValues(has_data, [0])
    error_string = "Data not found for the following required runsheet columns: "
    if len(listOfKeys) != 0:
        for key in listOfKeys:
            error_string = error_string + key + ', '
        raise ValueError(error_string)
    return

def getKeysByValues(dictOfElements, listOfValues):
    listOfKeys = list()
    listOfItems = dictOfElements.items()
    for item  in listOfItems:
        if item[1] in listOfValues:
            listOfKeys.append(item[0])
    return  listOfKeys 


def parse_range_list(rl):
    def collapse_range(ranges):
        end = None
        for value in ranges:
            yield range(max(end, value.start), max(value.stop, end)) if end else value
            end = max(end, value.stop) if end else value.stop
    def split_range(value):
        value = value.split(':')
        for val, prev in zip(value, chain((None,), value)):
            if val != '':
                val = int(val)
                if prev == '':
                    val *= -1
                yield val
    def parse_range(r):
        parts = list(split_range(r.strip()))
        if len(parts) == 0:
            return range(0, 0)
        elif len(parts) > 2:
            raise ValueError("Invalid range: {}".format(r))
        return range(parts[0], parts[-1] + 1)
    ranges = sorted(set(map(parse_range, rl.split(","))), key=lambda x: (x.start, x.stop))
    return chain.from_iterable(collapse_range(ranges))


# if __name__ == '__main__':
#   parser = argparse.ArgumentParser('A wrapper for running henipipe')
#   parser.add_argument('job', type=str, choices=['MAKERUNSHEET', 'ALIGN', 'NORM', 'SEACR'], help='a required string denoting segment of pipeline to run.  1) "MAKERUNSHEET" - to parse a folder of fastqs; 2) "ALIGN" - to perform alignment using bowtie and output bed files; 3) "NORM" - to normalize data to reference (spike in); 4) "SEACR" - to perform SEACR.')
#   parser.add_argument('--sample_flag', '-sf', type=str, help='FOR MAKERUNSHEET only string to identify samples of interest in a fastq folder')
#   parser.add_argument('--fastq_folder', '-fq', type=str, help='For MAKERUNSHEET only: Pathname of fastq folder (files must be organized in folders named by sample)')
#   parser.add_argument('--filter_high', '-fh', type=int, default=None, help='For ALIGN only: upper limit of fragment size to exclude, defaults is no upper limit.  OPTIONAL')
#   parser.add_argument('--filter_low', '-fl', type=int, default=None, help='For ALIGN only: lower limit of fragment size to exclude, defaults is no lower limit.  OPTIONAL')
#   parser.add_argument('--output', '-o', type=str, default=".", help='For MAKERUNSHEET only: Pathname to write runsheet.csv file (folder must exist already!!), Defaults to current directory')
#   parser.add_argument('--runsheet', '-r', type=str, help='tab-delim file with sample fields as defined in the script. - REQUIRED for all jobs except MAKERUNSHEET')
#   parser.add_argument('--log_prefix', '-l', type=str, default='henipipe.log', help='Prefix specifying log files for henipipe output from henipipe calls. OPTIONAL')
#   parser.add_argument('--select', '-s', type=str, default=None, help='To only run the selected row in the runsheet, OPTIONAL')
#   parser.add_argument('--debug', '-d', action='store_true', help='To print commands (For testing flow). OPTIONAL')
#   parser.add_argument('--bowtie_flags', '-b', type=str, default='--end-to-end --very-sensitive --no-mixed --no-discordant -q --phred33 -I 10 -X 700', help='For ALIGN: bowtie flags, OPTIONAL')
#   parser.add_argument('--cluster', '-c', type=str, default='PBS', choices=['PBS', 'SLURM'], help='Cluster software.  OPTIONAL Currently supported: PBS and SLURM')
#   parser.add_argument('--norm_method', '-n', type=str, default='coverage', choices=['coverage', 'read_count', 'spike_in'], help='For ALIGN and NORM: Normalization method, by "read_count", "coverage", or "spike_in".  If method is "spike_in", HeniPipe will align to the spike_in reference genome provided in runsheet. OPTIONAL')
#   parser.add_argument('--user', '-u', type=str, default=None, help='user for submitting jobs - defaults to username.  OPTIONAL')
#   parser.add_argument('--SEACR_norm', '-Sn', type=str, default='non', choices=['non', 'norm'], help='For SEACR: Normalization method; default is "non"-normalized, select "norm" to normalize using SEACR. OPTIONAL')
#   parser.add_argument('--SEACR_stringency', '-Ss', type=str, default='stringent', choices=['stringent', 'relaxed'], help='FOR SEACR: Default will run as "stringent", other option is "relaxed". OPTIONAL')
#   parser.add_argument('--verbose', '-v', default=False, action='store_true', help='Run with some additional ouput - not much though... OPTIONAL')
#   #call = '/home/sfurla/Scripts/runHeniPipe.py MAKERUNSHEET -sf mini -fq /active/furlan_s/Data/CNR/190801_CNRNotch/fastq/mini/fastqs'

#   #args = parser.parse_args(call.split(" ")[1:])
#   args = parser.parse_args()

#   #log
#   if args.debug == False:
#       LOGGER.info("Logging to %s... examine this file if samples fail." % args.log_prefix)

#   #deal with user
#   if args.user is None:
#       args.user = getpass.getuser()

#   #deal with paths
#   if args.job=="MAKERUNSHEET":
#       if os.path.isabs(args.fastq_folder) is False:
#           if args.fastq_folder == ".":
#               args.fastq_folder = os.getcwd()
#           else :
#               args.fastq_folder = os.path.abspath(args.fastq_folder)
#       if os.path.exists(args.fastq_folder) is False:
#           raise ValueError('Path: '+args.fastq_folder+' not found')
#       if os.path.isabs(args.output) is False:
#           if args.output == ".":
#               args.output = os.getcwd()
#           else :
#               args.output = os.path.abspath(args.output)
#       if os.path.exists(args.output) is False:
#           raise ValueError('Path: '+args.output+' not found')
#   if args.job != "MAKERUNSHEET":
#       if os.path.exists(args.runsheet) is False:
#           raise ValueError('Path: '+args.runsheet+' not found')

#   if args.job=="MAKERUNSHEET":
#       LOGGER.info("Parsing fastq folder - "+args.fastq_folder+" ...")
#       make_runsheet(folder=args.fastq_folder, output=args.output, furlan=True, sample_flag = args.sample_flag)
#       exit()

#   #parse and chech runsheet
#   args.runsheet = os.path.abspath(args.runsheet)
#   parsed_runsheet = list(parse_runsheet(args.runsheet))
#   check_runsheet(args, parsed_runsheet, verbose=args.verbose)

#   #deal with sample selection
#   if args.select is not None:
#       pare_down = [int(args.select) -1]
#   else:
#       pare_down = list(range(len(parsed_runsheet)))

#   if args.job=="ALIGN":
#       #deal with filtering
#       LOGGER.info("Aligning reads...")
#       Alignjob = Align(runsheet_data = [parsed_runsheet[i] for i in pare_down], debug=args.debug, cluster=args.cluster, bowtie_flags=args.bowtie_flags, log=args.log_prefix, user=args.user, norm_method=args.norm_method, filter = [args.filter_low, args.filter_high])
#       LOGGER.info("Submitting alignment jobs... Debug mode is %s" % args.debug)
#       Alignjob.run_job()

#   if args.job=="NORM":
#       LOGGER.info("Calculating %s", args.norm_method)
#       Normjob = Norm(runsheet_data = [parsed_runsheet[i] for i in pare_down], debug=args.debug, cluster=args.cluster, log=args.log_prefix, norm_method=args.norm_method, user=args.user)
#       LOGGER.info("Submitting bedgraph jobs... Debug mode is %s" % args.debug)
#       Normjob.run_job()

#   if args.job=="SEACR":
#       LOGGER.info("Running SEACR using settings: SEACR_norm = %s, SEACR_stringency = %s" % (args.SEACR_norm, args.SEACR_stringency))
#       SEACRjob = SEACR(runsheet_data = parsed_runsheet, pare_down = pare_down, debug=args.debug, cluster=args.cluster, norm=args.SEACR_norm, stringency=args.SEACR_stringency, user=args.user, log=args.log_prefix)
#       SEACRjob.run_job()


