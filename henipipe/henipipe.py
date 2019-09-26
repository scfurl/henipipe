#!/usr/bin/python
# PBS/SLURM cluster job submission in Python
#v 0.3.2

# A wrapper that performs the necessary pipeline for generation of CNR/CNT data
# Written by Scott Furlan with code inspiration from Andrew Hill's cellwrapper; 
# uses a custom python script samTobed.py which takes code from a fantastic 
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
import json
#import pandas as pd
#_ROOT = os.getcwd()
_ROOT = os.path.abspath(os.path.dirname(__file__))
GENOMES_JSON = os.path.join(_ROOT, 'data', 'genomes.json')
SEACR_SCRIPT = os.path.join(_ROOT, 'scripts', 'SEACR_1.1.sh')


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
                to_append = "#!/bin/bash\n#SBATCH --job-name=%s\n#SBATCH --output=outtmp\n#SBATCH --error=errtmp\n#SBATCH --ntasks=1\n%s\n{%s} 2>&1 | tee %s\nsed -e 's/^/[HENIPIPE] JOB: %s:\t\t/' %s >> %s\nrm %s\n" % (job_name, self.processor_line, command, log_file, job_name, log_file, self.log_name, log_file)
                #to_append = '#!/bin/bash\n#SBATCH --job-name=%s\n#SBATCH --ntasks=1\n%s\n%s' % (job_name, self.processor_line, command)
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
                    if (sys.version_info > (3, 0)):
                        proc.stdin.write(script.encode('utf-8'))
                        out, err = proc.communicate()
                    else:
                        proc.stdin.write(script)
                        out, err = proc.communicate()
                print(script)
                if self.debug==False:
                    print(out)
                    time.sleep(0.1)


class Align(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(Align, self).__init__(*args, **kwargs)
        self.bowtie_flags = kwargs.get('bowtie_flags')
        self.job = "HENIPIPE_ALIGN"
        self.pipe = not kwargs.get('no_pipe')
        self.threads = int(kwargs.get('threads'))
        self.gb_ram = int(kwargs.get('gb_ram'))
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
            sam2bed_string = """| samTobed - -o %s %s""" % (sample['bed_out']+'tmp', self.filter_string)
            if self.cluster=="SLURM":
                modules = """\nml bowtie2\nmodule load samtools\nmodule load Python/3.6.7-foss-2016b-fh1\n"""
            else:
                modules = """\nmodule load python\nmodule load bowtie2\nmodule load samtools\necho '\nRunning Bowtie piped to samTobed...\n[BOWTIE] Output:\n'\n"""
            norm_bowtie_flags='--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant -q --phred33 -I 10 -X 700'
            if self.pipe:
                commandline = """bowtie2 %s -p %s -1 %s -2 %s -x %s %s\n""" % (self.bowtie_flags, self.threads, fastq1, fastq2, sample['fasta'], sam2bed_string)
                commandline = commandline + """\necho 'Sorting Bed...\n'\nsort -k1,1 -k2n,2n %s > %s\n""" % (sample['bed_out']+'tmp', sample['bed_out'])
                commandline = commandline + """rm %s \n""" % (sample['bed_out']+'tmp')
            else:
                commandline = """bowtie2 %s -p %S -1 %s -2 %s -x %s -S %s\n""" % (self.bowtie_flags, self.threads, fastq1, fastq2, sample['fasta'], sample['sample']+".sam")
                commandline = commandline + """samtools view -bS %s > %s""" % (sample['sample']+".sam", sample['sample']+".bam")
                commandline = commandline + """\necho 'samToBed...\n'\nsamTobed %s -o %s %s""" % (sample['sample']+".sam", sample['bed_out']+'tmp', self.filter_string)
                commandline = commandline + """\necho 'Sorting Bed...\n'\nsort -k1,1 -k2n,2n %s > %s\n""" % (sample['bed_out']+'tmp', sample['bed_out'])
            if self.norm_method == "spike_in":
                commandline = commandline + """echo '\n[BOWTIE] Running Bowtie piped to samTobed.py for spikein... Output:\n'\nbowtie2 %s -p 4 -1 %s -2 %s -x %s | samTobed - -o %s\n""" % (norm_bowtie_flags, fastq1, fastq2, sample['spikein_fasta'], sample['spikein_bed_out']+'tmp')
                commandline = commandline + """\necho 'Sorting Bed for spikein...\n'sort -k1,1 -k2n,2n %s > %s\n""" % (sample['spikein_bed_out']+'tmp', sample['spikein_bed_out'])
                commandline = commandline + """rm %s \n""" % (sample['spikein_bed_out']+'tmp')
            commandline = modules + commandline
            command.append(commandline)
        return command


    def align_processor_line(self):
        if self.cluster=="PBS":
            return """select=1:mem=%sGB:ncpus=%s""" %(self.gb_ram*self.gb_ram, self.threads)
        if self.cluster=="SLURM":
            return '#SBATCH --cpus-per-task=%s\n#SBATCH --mem-per-cpu=%s000' %(self.threads, self.gb_ram)


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
                modules = """\nml bedtools\n"""
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
                modules = """\nsource /app/Lmod/lmod/lmod/init/bash\n"""
            else:
                modules = """\nmodule load R\nmodule load bedtools\n"""
            #commandline = """echo '\n[SEACR] Running SEACR... Output:\n'bash /home/sfurla/develop/SEACR/SEACR_1.1.sh %s %s %s %s %s""" % (sample['SEACR_in'], sample['SEACR_control'], self.norm, self.method, sample['SEACR_out'])
            commandline = """echo '\n[SEACR] Running SEACR... Output:\n'\nbash %s %s %s %s %s %s\n""" % (SEACR_SCRIPT, sample['SEACR_in'], sample['SEACR_control'], self.norm, self.method, sample['SEACR_out'])
            commandline = modules + commandline
            command.append(commandline)
        return command


    def SEACR_processor_line(self):
        if self.cluster=="PBS":
            return """select=1:mem=8GB:ncpus=2"""
        if self.cluster=="SLURM":
            return ''


class Merge(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(Merge, self).__init__(*args, **kwargs)
        self.job = "HENIPIPE_MERGE"
        self.out = kwargs.get('out')
        self.runsheet_data = self.Merge_match(pare_down = kwargs.get('pare_down'))
        self.processor_line = self.Merge_processor_line()
        self.command = self.Merge_executable(pare_down = kwargs.get('pare_down'))
        self.script = self.generate_job()
    def __call__():
        pass

    def Merge_match(self, pare_down):
        key_data = [self.runsheet_data[i].get("merge_key") for i in pare_down]
        bg_data = [self.runsheet_data[i].get("bedgraph") for i in pare_down]
        merge_dict = dict.fromkeys(key_data, "NotFound")
        samples = []
        for key in merge_dict.keys():
            # do something with value
            samples.append({    "sample" : key,
                                "files_to_merge": list(compress(bg_data, is_in(key, key_data))),})
        return(samples)


    def Merge_executable(self, pare_down):
        commandline=""
        command = []
        #print("Runmode is " + self.runmode)
        #print(keys)
        for i in self.runsheet_data:
            #print(key)
            seperator = ' '
            nfiles = len(i.get("files_to_merge"))
            bedgraph_line = seperator.join(i.get("files_to_merge"))
            bedgraph_out=str(os.path.join(self.out, i.get("sample")))+"_merged.bedgraph"
            JOBSTRING = self.id_generator(size=10)
            if self.cluster=="SLURM":
                modules = """\nsource /app/Lmod/lmod/lmod/init/bash\nmodule load bedtools\n"""
            else:
                modules = """\nmodule load bedtools\n"""
            commandline = """echo '\n[MERGE] Merging bedgraphs:\n%s'\nbedtools unionbedg -i %s | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $0"\t"sum/(NF-4+1); }' > %s\nsleep 10\ncut -d$'\t' -f1-3,%s %s > %s\nrm %s\n""" % (bedgraph_line, bedgraph_line, bedgraph_out+'temp', (nfiles+4), bedgraph_out+'temp', bedgraph_out, bedgraph_out+'temp')
            commandline = modules + commandline
            command.append(commandline)
        return command


    def Merge_processor_line(self):
        if self.cluster=="PBS":
            return """select=1:mem=8GB:ncpus=2"""
        if self.cluster=="SLURM":
            return ''


class MACS2(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(MACS2, self).__init__(*args, **kwargs)
        self.job = "HENIPIPE_MACS2"
        self.merged = kwargs.get('merged')
        self.norm = kwargs.get('norm')
        self.runsheet_data = self.MACS2_match(pare_down = kwargs.get('pare_down'))
        print(self.runsheet_data)
        self.processor_line = self.MACS2_processor_line()
        self.command = self.MACS2_executable(pare_down = kwargs.get('pare_down'))
        self.script = self.generate_job()
    def __call__():
        pass

    def MACS2_match(self, pare_down):
        #will need to change this when multiple selections are implemented; for now just allow user to specify sample, then find control
        # if self.merged:
        #     desired_samples = [self.runsheet_data[i] for i in pare_down]
        #     #desired_samples = [parsed_runsheet[i] for i in pare_down]
        #     key_data = [i.get("merge_key") for i in desired_samples]
        #     match_data = [i.get("merge_MACS2_key") for i in desired_samples]
        #     unique_keys = unique(key_data)
        #     run_list = []
        #     for key in unique_keys:
        #         #find out if file is sample or control by searching lists 
        #         query = [match_data[i] for i in which(key, key_data)]
        #         bools = [bool(re.search(r'._CONTROL$', i)) for i in query]
        #         is_control = all_the_same(bools)
        #         if is_control == 'mixed':
        #             raise ValueError("Some discrepency between merge_key and MACS2_key ")
        #         if is_control:
        #             control_filename = key +"_merged.bedgraph"
        #             sample = re.sub("_CONTROL", "", query[0])
        #             sample_filename = key_data[which(sample, match_data)[0]]+"_merged.bedgraph"
        #             run_list.append({   "MACS2_in": sample_filename,
        #                                 "MACS2_control": control_filename,
        #                                 "sample": sample_filename})
        #     return(run_list)
        if self.merged:
            desired_samples = [self.runsheet_data[i] for i in pare_down]
            #desired_samples = [parsed_runsheet[i] for i in pare_down]
            key_data = [i.get("merge_key") for i in desired_samples]
            match_data = [i.get("SEACR_key") for i in desired_samples]
            unique_keys = unique(key_data)
            run_list = []
            for key in unique_keys:
                #find out if file is sample or control by searching lists 
                query = [match_data[i] for i in which(key, key_data)]
                bools = [bool(re.search(r'._CONTROL$', i)) for i in query]
                is_control = all_the_same(bools)
                if is_control == 'mixed':
                    raise ValueError("Some discrepency between merge_key and SEACR_key ")
                if is_control:
                    control_filename = key +"_merged.bedgraph"
                    sample = re.sub("_CONTROL", "", query[0])
                    sample_filename = key_data[which(sample, match_data)[0]]+"_merged.bedgraph"
                    run_list.append({   "treatment_in": sample_filename,
                                        "control_in": control_filename,
                                        "sample": sample_filename})
            return(run_list)
        else:
            desired_samples = [self.runsheet_data[i] for i in pare_down]
            sk = [i.get('MACS2_key') for i in desired_samples]
            controls_b = [bool(re.search(r'._CONTROL$', i)) for i in sk]
            controls = list(compress(desired_samples, controls_b))
            samples_b = [not i for i in controls_b]
            samples = list(compress(desired_samples, samples_b))
            for sample in samples:
                control_name = sample.get('MACS2_key')+"_CONTROL"
                control_bed = next(item for item in controls if item["MACS2_key"] == control_name).get('bedgraph')
                sample.update( {'MACS2_in' : sample.get('bedgraph')})
                sample.update( {'MACS2_control' : control_bed})
            return samples



    def MACS2_executable(self, pare_down):
        commandline=""
        command = []
        #print("Runmode is " + self.runmode)
        for item in self.runsheet_data:
            JOBSTRING = self.id_generator(size=10)
            #print(sample)
            macs2_out = re.sub("_merged.bedgraph", "_MACS2.bedgraph",item['MACS2_in'])
            if self.cluster=="SLURM":
                modules = """\nsource /app/Lmod/lmod/lmod/init/bash\nmodule load MACS2"""
            else:
                modules = """\n"""
            #commandline = """echo '\n[SEACR] Running SEACR... Output:\n'bash /home/sfurla/develop/SEACR/SEACR_1.1.sh %s %s %s %s %s""" % (sample['SEACR_in'], sample['SEACR_control'], self.norm, self.method, sample['SEACR_out'])
            commandline = """echo '\n[MACS2] Running MACS2... Output:\n'\nmacs2 bdgcmp -t %s -c %s -o %s -m FE\n""" % (item['MACS2_in'], item['MACS2_control'], macs2_out)
            commandline = modules + commandline
            command.append(commandline)
        return command

    def MACS2_processor_line(self):
        if self.cluster=="PBS":
            return """select=1:mem=8GB:ncpus=2"""
        if self.cluster=="SLURM":
            return ''

def which(string, list):
    return [i for i, x in enumerate(list) if x == string]

def unique_indices(list):
    result = {}
    for key, val in zip(keys, vals):
        if key not in result:
            result[key] = 0
        result[key] += val
    return(result)


def unique(list1):
    # function to get unique values
    # intilize a null list
    unique_list = []
    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    return(unique_list)

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

def is_in(string, list):
    outlist=[]
    #This function returns a list of booleans for matching a string in a list of strings
    for i in list:
        if i == string:
            outlist.append(True)
        else:
            outlist.append(False)
    return(outlist)

def reduce_concat(x, sep=""):
    return functools.reduce(lambda x, y: str(x) + sep + str(y), x)

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

def all_the_same(list_in):
    #this function takes a list of booleans and characterizes them as all true, all false or mixed.
    if all(list_in) is True:
        return(True)
    elif all([not i for i in list_in]) is True:
        return(False)
    else:
        return("mixed")

def find_colnames(runsheet, header=True):
    if header==True:
        with open(runsheet, 'r') as f:
            reader = csv.reader(f)
            return(next(reader))            # read header


def load_genomes(genomes_file):
    with open(genomes_file, "r") as read_file:
        genome_data = json.load(read_file)
    return genome_data

def make_runsheet(folder, sample_flag, genome_key, output="./henipipeout", fasta=None, spikein_fasta=None, genome_sizes=None):
    genome_data = load_genomes(GENOMES_JSON).get(genome_key)
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
            'MERGE_key': i.get('directory_short'), \
            'SEACR_out': os.path.join(output, i.get('directory_short')+"_SEACR"), \
            'fasta': genome_data.get('fasta'), 'spikein_fasta': genome_data.get('spikein_fasta'), 'genome_sizes':  genome_data.get('genome_sizes')})
    #print(good_dat)
    #keys = good_dat[0].keys()
    keys = ["sample", "SEACR_key", "MERGE_key", "fasta", "spikein_fasta", "genome_sizes", "fastq1", "fastq2", "bed_out", "spikein_bed_out", "bedgraph", "SEACR_out"]
    with open(os.path.join(output, 'runsheet.csv'), 'w') as output_file:
        dict_writer = csv.DictWriter(output_file, fieldnames = keys, extrasaction='ignore')
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
            if verbose: print("header present, but no or incomplete data for "+parameter)
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
