import sys
import argparse
import logging
import getpass
import os
from . import samTobed
from . import pyWriter
from . import henipipe

POLL_TIME = 5
LOG_PREFIX = '[HENIPIPE]: '

# Set up a basic logger
LOGGER = logging.getLogger('something')
myFormatter = logging.Formatter('%(asctime)s: %(message)s')
handler = logging.StreamHandler()
handler.setFormatter(myFormatter)
LOGGER.addHandler(handler)
LOGGER.setLevel(logging.DEBUG)
myFormatter._fmt = "[HENIPIPE]: " + myFormatter._fmt


def run_henipipe(args=None):
    if args is None:
        args = sys.argv[1:]
    parser = argparse.ArgumentParser('A wrapper for running henipipe')
    parser.add_argument('job', type=str, choices=['SC', 'MAKERUNSHEET', 'ALIGN', 'SCALE', 'MERGE', 'SEACR', 'MACS2', 'AUC', 'GENOMESFILE', 'FASTQC', 'TRIM', 'BIGWIG'], help='a required string denoting segment of pipeline to run.  1) "MAKERUNSHEET" - to parse a folder of fastqs; 2) "ALIGN" - to perform alignment using bowtie and output bed files; 3) "SCALE" - to normalize data to reference (spike in); 4) "MERGE" - to merge bedgraphs 5) "SEACR" - to perform SEACR; 6) "MACS" - to perform MACS2; 7) "AUC" - to calculate AUC between normalized bedgraph using a peak file; 8) "GENOMESFILE" - print location of genomes.json file; 9) "FASTQC" - run fastqc on cluster; 10) run trimmotatic on cluster; 11) make Bigwigs from bedgraphs;')
    parser.add_argument('--sample_flag', '-sf', type=str, default="", help='FOR MAKERUNSHEET only string to identify samples of interest in a fastq folder')
    parser.add_argument('--fastq_folder', '-fq', type=str, help='For SC and MAKERUNSHEET only: Pathname of fastq folder (files must be organized in folders named by sample)')
    parser.add_argument('--trim_folder', '-tf', type=str, default = ".", help='REQURIED, For TRIM only: Pathname of output folder; Note that all trimmed fastqs will be placed in the same folder')
    parser.add_argument('--trim_args', '-ta', type=str, default = "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36", help='OPTIONAL, For TRIM only: Args to pass to trimmomatic')
    parser.add_argument('--organized_by', '-by', type=str, choices=['folder', 'file'], default='folder', help='Option to specify how fastq or unbam folder is organized')
    parser.add_argument('--genome_key', '-gk', type=str, default="default", help='For MAKERUNSHEET only: abbreviation to use "installed" genomes in the runsheet (See README.md for more details')
    parser.add_argument('--split_char', '-sc', type=str, default="_R1_", help='Character by which to split the fastqfile name into samples, OPTIONAL and for MAKERUNSHEET only')
    parser.add_argument('--R1_char', '-r1c', type=str, default="_R1_", help='Character by which to split the fastqfile name into read1, OPTIONAL and for MAKERUNSHEET only; default = "_R1_"')
    parser.add_argument('--R2_char', '-r2c', type=str, default="_R2_", help='Character by which to split the fastqfile name into read2, OPTIONAL and for MAKERUNSHEET only; default = "_R2_"')
    parser.add_argument('--ext', '-e', type=str, default=".fastq.gz", help='suffix of fastq files, OPTIONAL and for MAKERUNSHEET only')
    parser.add_argument('--filter_high', '-fh', type=int, default=None, help='For ALIGN only: upper limit of fragment size to exclude, defaults is no upper limit.  OPTIONAL')
    parser.add_argument('--filter_low', '-fl', type=int, default=None, help='For ALIGN only: lower limit of fragment size to exclude, defaults is no lower limit.  OPTIONAL')
    parser.add_argument('--output', '-o', type=str, default=".", help='Pathname to output folder (note this folder must exist already!!), Defaults to current directory')
    parser.add_argument('--runsheet', '-r', type=str, default="runsheet.csv", help='tab-delim file with sample fields as defined in the script. - REQUIRED for all jobs except MAKERUNSHEET')
    parser.add_argument('--log_prefix', '-l', type=str, default='henipipe.log', help='Prefix specifying log files for henipipe output from henipipe calls. OPTIONAL')
    parser.add_argument('--select', '-s', type=str, default=None, help='To only run the selected row in the runsheet, OPTIONAL')
    parser.add_argument('--debug', '-d', action='store_true', help='To print commands (For testing flow). OPTIONAL')
    parser.add_argument('--bowtie_flags', '-b', type=str, default='--end-to-end --very-sensitive --no-mixed --no-discordant -q --phred33 -I 10 -X 700', help='For ALIGN: bowtie flags, OPTIONAL')
    parser.add_argument('--cluster', '-c', type=str, default='SLURM', choices=['PBS', 'SLURM', 'local'], help='Cluster software.  OPTIONAL Currently supported: PBS, SLURM and local')
    parser.add_argument('--threads', '-t', type=str, default="8", help='number of threads; default: 8')
    parser.add_argument('--gb_ram', '-gb', type=str, default=None, help='gigabytes of RAM per thread')
    parser.add_argument('--install', '-i', type=str, default=None, help='FOR GENOMESFILE: location of file to install as a new genomes.json file, existing genomes.json will be erased')
    parser.add_argument('--norm_method', '-n', type=str, default='coverage', choices=['coverage', 'read_count', 'spike_in', 'none'], help='For ALIGN and SCALE: Normalization method, by "read_count", "coverage", or "spike_in".  If method is "spike_in", HeniPipe will align to the spike_in reference genome provided in runsheet. OPTIONAL')
    parser.add_argument('--user', '-u', type=str, default=None, help='user for submitting jobs - defaults to username.  OPTIONAL')
    parser.add_argument('--SEACR_norm', '-Sn', type=str, default='non', choices=['non', 'norm'], help='For SEACR: Normalization method; default is "non"-normalized, select "norm" to normalize using SEACR. OPTIONAL')
    parser.add_argument('--SEACR_fdr', '-Sf', type=str, default='0.05',  help='For SEACR: Used to set FDR threshold when control is not used. OPTIONAL')
    parser.add_argument('--SEACR_stringency', '-Ss', type=str, default='stringent', choices=['stringent', 'relaxed'], help='FOR SEACR: Default will run as "stringent", other option is "relaxed". OPTIONAL')
    parser.add_argument('--keep_files', '-k', action ='store_true', default=False, help='FOR ALIGN: use this flag to turn off piping (Will generate all files).')
    parser.add_argument('--verbose', '-v', default=False, action='store_true', help='Run with some additional ouput - not much though... OPTIONAL')
    """
    call = 'henipipe MAKERUNSHEET -fq ../fastq -sf mini -gk heni_hg38 -o .'
    call = 'henipipe MACS2 -r ./runsheet.csv -d -mk -s 1:10'
    call = 'henipipe GENOMESFILE'
    call = 'henipipe MAKERUNSHEET -fq ../fastq'
    call = 'henipipe MAKERUNSHEET -fq ../fastq'
    call = 'henipipe ALIGN -r runsheet.csv -d'
    args = parser.parse_args(call.split(" ")[1:])
    """
    args = parser.parse_args()

    #deal with user
    if args.user is None:
        args.user = getpass.getuser()


    if args.job=="SC":
        if os.path.isabs(args.fastq_folder) is False:
            if args.fastq_folder == ".":
                args.fastq_folder = os.getcwd()
            else :
                args.fastq_folder = os.path.abspath(args.fastq_folder)
        if os.path.exists(args.fastq_folder) is False:
            raise ValueError('Path: '+args.fastq_folder+' not found')
        if os.path.isabs(args.output) is False:
            if args.output == ".":
                args.output = os.getcwd()
            else :
                args.output = os.path.abspath(args.output)
        if os.path.exists(args.output) is False:
            raise ValueError('Path: '+args.output+' not found')
        LOGGER.info("Running Single Cell processing...")
        if args.select is not None:
            select = [i for i in list(henipipe.parse_range_list(args.select))]
        else:
            select = None
        SCjob = henipipe.SCAlign(folder=args.fastq_folder, output=args.output,
            threads = args.threads, ram = args.gb_ram, debug=args.debug, sample_flag = args.sample_flag,
            genome_key = args.genome_key, no_pipe=args.keep_files, cluster=args.cluster,
            bowtie_flags=args.bowtie_flags, log=args.log_prefix, user=args.user,
            ext=args.ext,  r1_char=args.R1_char, strsplit= args.split_char,
            r2_char=args.R2_char, fname=args.runsheet, organized_by=args.organized_by,
            filter = [args.filter_low, args.filter_high], select = select)
        dependency = SCjob.run_dependency_job()
        SCjob = henipipe.SCMerge(runsheet_data = SCjob.runsheet_data, output=args.output, debug=args.debug,
            dependency = dependency, cluster=SCjob.environs.cluster, log=SCjob.environs.log, user=SCjob.environs.user,
            genome_key = args.genome_key, norm=args.SEACR_norm, stringency=args.SEACR_stringency, fdr_thresh = args.SEACR_fdr)
        SCjob.run_dependency_job()
        exit()

    #log


    if args.job=="GENOMESFILE":
        _ROOT = os.path.abspath(os.path.dirname(__file__))
        if args.install is None:
            GENOMES_JSON = os.path.join(_ROOT, 'data', 'genomes.json')
            print(GENOMES_JSON)
        if args.install is not None:
            from shutil import copyfile
            args.install = os.path.abspath(args.install)
            copyfile(args.install, os.path.join(_ROOT, 'data', 'genomes.json'))
        exit()
    #log


    #deal with paths
    if args.job=="MAKERUNSHEET":
        if os.path.isabs(args.fastq_folder) is False:
            if args.fastq_folder == ".":
                args.fastq_folder = os.getcwd()
            else :
                args.fastq_folder = os.path.abspath(args.fastq_folder)
        if os.path.exists(args.fastq_folder) is False:
            raise ValueError('Path: '+args.fastq_folder+' not found')
        if os.path.isabs(args.output) is False:
            if args.output == ".":
                args.output = os.getcwd()
            else :
                args.output = os.path.abspath(args.output)
        if os.path.exists(args.output) is False:
            raise ValueError('Path: '+args.output+' not found')
    if args.job != "MAKERUNSHEET":
        if os.path.exists(args.runsheet) is False:
            raise ValueError('Path: '+args.runsheet+' not found')
        args.output = os.path.abspath(args.output)


    if args.job=="MAKERUNSHEET":
        LOGGER.info("Parsing fastq folder - "+args.fastq_folder+" ...")
        LOGGER.info("Writing runsheet to - "+os.path.join(args.output, 'runsheet.csv')+" ...")
        LOGGER.info("Using genome_key - "+args.genome_key+" ...")
        #henipipe.make_runsheet(folder=args.fastq_folder, output=args.output, sample_flag = args.sample_flag, genome_key = args.genome_key, no_pipe=args.keep_files)
        henipipe.make_runsheet(folder=args.fastq_folder, output=args.output, sample_flag = args.sample_flag, genome_key = args.genome_key, \
            no_pipe=args.keep_files, ext=args.ext,  r1_char=args.R1_char, strsplit= args.split_char, \
            r2_char=args.R2_char, fname=args.runsheet, organized_by=args.organized_by)
        exit()

    #parse and chech runsheet
    args.runsheet = os.path.abspath(args.runsheet)

    """
    parsed_runsheet = list(parse_runsheet(args.runsheet))
    check_runsheet(args, parsed_runsheet, verbose=args.verbose)
    """
    parsed_runsheet = list(henipipe.parse_runsheet(args.runsheet))

    henipipe.check_runsheet(args, parsed_runsheet, verbose=args.verbose)

    #deal with sample selection
    if args.select is not None:
        parsed_runsheet = [parsed_runsheet[i-1] for i in list(henipipe.parse_range_list(args.select))]

    if args.debug == False:
        LOGGER.info("Logging to %s... examine this file if samples fail." % args.log_prefix)

    if args.job=="FASTQC":
        LOGGER.info("Running fastqc on all fastqs in runsheet")
        Fastqcjob = henipipe.Fastqc(runsheet_data = parsed_runsheet, threads = args.threads, gb_ram = args.gb_ram, debug=args.debug, cluster=args.cluster, log=args.log_prefix, user=args.user)
        Fastqcjob.run_job()
        exit()

    if args.job=="TRIM":
        if os.path.isabs(args.trim_folder) is False:
            if args.trim_folder == ".":
                args.trim_folder = os.getcwd()
            else :
                args.trim_folder = os.path.abspath(args.trim_folder)
        if os.path.exists(args.trim_folder) is False:
            raise ValueError('Path: '+args.trim_folder+' not found')
        LOGGER.info("Running trimmotatic on all fastqs in runsheet")
        Trimjob = henipipe.Trim(runsheet_data = parsed_runsheet, trim_args = args.trim_args, trim_folder  = args.trim_folder, threads = args.threads, gb_ram = args.gb_ram, debug=args.debug, cluster=args.cluster, log=args.log_prefix, user=args.user)
        Trimjob.run_job()
        exit()

    if args.job=="ALIGN":
        #deal with filtering
        LOGGER.info("Aligning reads...")




        #Alignjob = Align(runsheet_data = parsed_runsheet, threads = args.threads, gb_ram = args.gb_ram, debug=args.debug, no_pipe=args.keep_files, cluster=args.cluster, bowtie_flags=args.bowtie_flags, log=args.log_prefix, user=args.user, norm_method=args.norm_method, filter = [args.filter_low, args.filter_high])


        Alignjob = henipipe.Align(runsheet_data = parsed_runsheet, threads = args.threads, gb_ram = args.gb_ram, debug=args.debug, no_pipe=args.keep_files, cluster=args.cluster, bowtie_flags=args.bowtie_flags, log=args.log_prefix, user=args.user, norm_method=args.norm_method, filter = [args.filter_low, args.filter_high])
        LOGGER.info("Submitting alignment jobs... Debug mode is %s" % args.debug)
        Alignjob.run_job()
        exit()

    if args.job=="SCALE":
        LOGGER.info("Calculating %s", args.norm_method)
        Scalejob = henipipe.Scale(runsheet_data = parsed_runsheet, threads = args.threads, gb_ram = args.gb_ram, debug=args.debug, cluster=args.cluster, log=args.log_prefix, norm_method=args.norm_method, user=args.user)
        LOGGER.info("Submitting bedgraph jobs... Debug mode is %s" % args.debug)
        Scalejob.run_job()
        exit()

    if args.job=="BIGWIG":
        LOGGER.info("Making Bigwigs:")
        Bigwigjob = henipipe.Bigwig(runsheet_data = parsed_runsheet, debug=args.debug, cluster=args.cluster, log=args.log_prefix, norm_method=args.norm_method, user=args.user)
        LOGGER.info("Submitting bigwig jobs... Debug mode is %s" % args.debug)
        Bigwigjob.run_job()
        exit()

    if args.job=="MERGE":
        Mergejob = henipipe.Merge(runsheet_data = parsed_runsheet, threads = args.threads, gb_ram = args.gb_ram, debug=args.debug, cluster=args.cluster, log=args.log_prefix, norm_method=args.norm_method, user=args.user, out=args.output)
        #Mergejob = Merge(runsheet_data = parsed_runsheet, debug=args.debug, cluster=args.cluster, log=args.log_prefix, norm_method=args.norm_method, user=args.user)
        LOGGER.info("Submitting merge-bedgraph jobs... Debug mode is %s" % args.debug)
        Mergejob.run_job()
        exit()

    if args.job=="SEACR":
        LOGGER.info("Running SEACR using settings: SEACR_norm = %s, SEACR_stringency = %s" % (args.SEACR_norm, args.SEACR_stringency))
        SEACRjob = henipipe.SEACR(runsheet_data = parsed_runsheet, threads = args.threads, gb_ram = args.gb_ram, debug=args.debug, cluster=args.cluster, norm=args.SEACR_norm, stringency=args.SEACR_stringency, user=args.user, log=args.log_prefix)
        SEACRjob.run_job()
        exit()

    if args.job=="MACS2":
        LOGGER.info("Running MACS2")
        MACS2job = henipipe.MACS2(runsheet_data = parsed_runsheet, threads = args.threads, gb_ram = args.gb_ram, debug=args.debug, cluster=args.cluster, user=args.user, log=args.log_prefix, out=args.output)
        MACS2job.run_job()
        exit()

    if args.job=="AUC":
        LOGGER.info("Running AUC")
        AUCjob = henipipe.AUC(runsheet_data = parsed_runsheet, threads = args.threads, gb_ram = args.gb_ram, debug=args.debug, no_pipe=args.keep_files, cluster=args.cluster, user=args.user, log=args.log_prefix, out=args.output, norm=args.SEACR_norm, stringency=args.SEACR_stringency)
        AUCjob.run_job()
        exit()


if __name__ == "__main__":
    run_henipipe()

"""
[parsed_runsheet[i-1] for i in list(parse_range_list("1:4,11,12"))]

"""



