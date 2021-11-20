#!/usr/bin/python
# PBS/SLURM cluster job submission in Python

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
import pathlib
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
from datetime import datetime


#_ROOT = os.getcwd()
_ROOT = os.path.abspath(os.path.dirname(__file__))
#_ROOT = '/Users/sfurla/Box Sync/PI_FurlanS/computation/develop/henipipe/henipipe'
GENOMES_JSON = os.path.join(_ROOT, 'data', 'genomes.json')
SEACR_SCRIPT = os.path.join(_ROOT, 'scripts', 'SEACR_1.3.sh')
ENVIRONS_JSON = os.path.join(_ROOT, 'data', 'environs.json')



class SampleFactory:
    def __init__(self, *args, **kwargs):
        self.environs = Environs(cluster = kwargs.get('cluster'), user = kwargs.get('user'), log = kwargs.get('log'), threads = kwargs.get('threads'), gb_ram = kwargs.get('gb_ram'))

        #remove later
        self.user = kwargs.get('user')
        self.cluster = kwargs.get('cluster')
        self.threads = kwargs.get('threads')
        self.gb_ram = kwargs.get('gb_ram')
        self.cluster = kwargs.get('cluster')
        self.runsheet_data = kwargs.get('runsheet_data')
        self.debug = kwargs.get('debug')
        self.dependency = kwargs.get('dependency')
        # self.log_name = kwargs.get('log')
    def __call__():
        pass

    def id_generator(self, size=10, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

    def run_job(self):
        popen_command = self.environs.popen_command
        with open('out.log', mode='w') as out_log, open('err.log', mode='w') as err_log:
            for script in self.bash_scripts:
                if self.debug==False:
                    # Open a pipe to the command.
                    proc = Popen(popen_command, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
                    if (sys.version_info > (3, 0)):
                        proc.stdin.write(script.encode('utf-8'))
                        out, err = proc.communicate()
                    else:
                        proc.stdin.write(script)
                        out, err = proc.communicate()
                # Print your job and the system response to the screen as it's submitted
                print(script)
                if self.debug==False:
                    #print(err)
                    print(str(out))
                    out_log.write(str(out))
                    err_log.write(str(err))
                time.sleep(0.1)

    def run_dependency_job(self):
        popen_command = self.environs.popen_command
        for script in self.bash_scripts:
            if self.debug==False:
                # Open a pipe to the command.
                proc = Popen(popen_command, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
                if (sys.version_info > (3, 0)):
                    proc.stdin.write(script.encode('utf-8'))
                    out, err = proc.communicate()
                else:
                    proc.stdin.write(script)
                    out, err = proc.communicate()
            # Print your job and the system response to the screen as it's submitted
            print(script)
            if self.debug==False:
                #print(err)
                print(out.decode("utf-8"))
                #out_log.write(str(out))
                #err_log.write(str(err))
                waiting = re.sub("Submitted batch job ", "", out.decode("utf-8")).strip()
                return(waiting)
            time.sleep(0.1)
            return("No jobs submitted")




class Environs:
    def __init__(self, *args, **kwargs):
        self.cluster = kwargs.get('cluster')
        self.user = kwargs.get('user')
        self.log = kwargs.get('log')
        self.environs_data = self.load_environs(ENVIRONS_JSON).get(self.cluster)
        self.popen_command = self.environs_data["popen"]
        self.threads = kwargs.get('threads')
        self.ram = kwargs.get('gb_ram')
        self.array_job=False
        self.dependency_job=None

    def id_generator(self, size=10, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

    def load_environs(self, environs_file):
        with open(environs_file, "r") as read_file:
            data = json.load(read_file)
        return data

    def get_processor_line(self, *args, **kwargs):
        if kwargs.get('threads') is None:
            threads = self.environs_data["resources"][kwargs.get('job')]["threads"]
        else:
            threads = kwargs.get('threads')
        if kwargs.get('gb_ram') is None:
            ram = self.environs_data["resources"][kwargs.get('job')]["ram"]
        else:
            ram = kwargs.get('ram')
        return [threads, ram]

    def generate_job(self, commands, job):
        bash_scripts=[]
        torun = len(commands)
        #print("Array count = "+str(self.array_count))
        threads,ram = self.get_processor_line(threads = self.threads, ram = self.ram, job = job)
        if self.array_job:
            bash_script = self.assemble_script(LOG_FILE = self.log,
                                    MODULES = self.environs_data["resources"][job]["modules"],
                                    TEMP_LOG_FILE = self.id_generator(),
                                    JOB_NAME = (job + "_ARRAY"),
                                    COMMAND = commands[0][1],
                                    ARRAY_COUNT = str(self.array_count),
                                    RAM = ram,
                                    NUM_JOBS = torun,
                                    THREADS = threads,
                                    USER = self.user)
                                    # TIME = str(datetime.now().strftime("%Y-%m-%d_%H:%M:%S")))
            bash_scripts.append(bash_script)
        else:
            for i in range(torun):
            # bash_script = {  "PBS" : "#!/bin/bash\n#PBS -N %s\n#PBS -l %s\n#PBS -j oe\n#PBS -o $PBS_O_WORKDIR/logtmp\n#PBS -A %s\ncd $PBS_O_WORKDIR\n{%s} 2>&1 | tee %s\nsed -e 's/^/[HENIPIPE] JOB: %s:\t\t/' %s >> %s\nrm %s\n" % (job_name, self.processor_line, self.user, command, log_file, job_name, log_file, self.log_name, log_file),
                            # "SLURM" : "#!/bin/bash\n#SBATCH --job-name=%s\n#SBATCH --output=tmp\n#SBATCH --error=tmp\n#SBATCH --ntasks=1\n%s\n{%s} 2>&1 | tee %s\nsed -e 's/^/[HENIPIPE] JOB: %s:\t\t/' %s >> %s\nrm %s\n" % (job_name, self.processor_line, command, log_file, job_name, log_file, self.log_name, log_file)}
                bash_script = self.assemble_script(LOG_FILE = self.log,
                                        MODULES = self.environs_data["resources"][job]["modules"],
                                        TEMP_LOG_FILE = self.id_generator(),
                                        JOB_NAME = (job + "_" + commands[i][0]),
                                        COMMAND = commands[i][1],
                                        RAM = ram,
                                        THREADS = threads,
                                        USER = self.user)
                                        # TIME = str(datetime.now().strftime("%Y-%m-%d_%H:%M:%S")))
                bash_scripts.append(bash_script)
        return bash_scripts

    def assemble_script(self, *args, **kwargs):
        script_list=[]
        order=[]
        fn_args = kwargs
        if self.array_job:
            which_header = "array_script_lines"
        else:
            which_header = "script_lines"
        for key, value in self.environs_data[which_header].items():
            script_list.append(value)
            order.append(key)
        script_list = [x for _,x in sorted(zip(order,script_list))]
        lines_unparsed = [x[0] for x in script_list]
        values_to_insert = [x[1].split("|") for x in script_list]
        lines_parsed = []
        if self.cluster == 'local':
            lines_parsed = [kwargs['COMMAND']]
        else:
            global to_test
            to_test=[lines_unparsed, values_to_insert, fn_args]
            for i in range(len(lines_unparsed)):
                if values_to_insert[i][0] is "":
                    lines_parsed.append(lines_unparsed[i])
                else:
                    string=lines_unparsed[i]
                    for j in range(len(values_to_insert[i])):
                        string = re.sub("<--{0}-->".format(j), fn_args.get(values_to_insert[i][j]), string)
                    lines_parsed.append(string)
        if self.dependency_job is not None:
            line = "#SBATCH --depend=afterok:%s" % (self.dependency_job)
            lines_parsed.insert(1, line)
        return "\n".join(lines_parsed)


class Fastqc(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(Fastqc, self).__init__(*args, **kwargs)
        self.job = "HENIPIPE_FASTQC"
        self.commands = self.fastqc_executable()
        self.bash_scripts = self.environs.generate_job(self.commands, self.job)
    def __call__():
        pass

    def fastqc_executable(self):
        commandline=""
        command = []
        for sample in self.runsheet_data:
            fastq1=re.sub('\t', ' ', sample['fastq1'])
            fastq2=re.sub('\t', ' ', sample['fastq2'])
            JOBSTRING = self.id_generator(size=10)
            commandline = """fastqc %s %s\n""" % (fastq1, fastq2)
            command.append([sample['sample'], commandline])
        return command


class Trim(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(Trim, self).__init__(*args, **kwargs)
        self.job = "HENIPIPE_TRIM"
        self.trim_args = kwargs.get('trim_args')
        self.trim_folder = kwargs.get('trim_folder')
        self.commands = self.trim_executable()
        self.bash_scripts = self.environs.generate_job(self.commands, self.job)
    def __call__():
        pass

    def trim_executable(self):
        commandline=""
        command = []
        for sample in self.runsheet_data:
            fastq1=sample['fastq1']
            fastq2=sample['fastq2']
            #make unpaired folder
            unpaired_folder = os.path.join(pathlib.Path(self.trim_folder), "unpaired")
            if os.path.exists(unpaired_folder) is False:
                os.mkdir(unpaired_folder)
            pairedoutfastq1=change_dirname(pathlib.Path(sample['fastq1']), self.trim_folder)
            pairedoutfastq2=change_dirname(pathlib.Path(sample['fastq2']), self.trim_folder)
            unpairedoutfastq1=change_dirname(pathlib.Path(sample['fastq1']), unpaired_folder)
            unpairedoutfastq2=change_dirname(pathlib.Path(sample['fastq2']), unpaired_folder)
            JOBSTRING = self.id_generator(size=10)
            commandline = """java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE %s %s %s %s %s %s %s\n""" % (fastq1, fastq2, pairedoutfastq1, pairedoutfastq2, unpairedoutfastq1, unpairedoutfastq2, self.trim_args)
            command.append([sample['sample'], commandline])
        return command

class SCAlign(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(SCAlign, self).__init__(*args, **kwargs)
        self.bowtie_flags = kwargs.get('bowtie_flags')
        self.job = "HENIPIPE_SCALIGN"
        self.pipe = not kwargs.get('no_pipe')
        self.filter_string = self.make_filter_string(kwargs.get('filter')[0], kwargs.get('filter')[1])
        self.genome_data = load_genomes(GENOMES_JSON).get(kwargs.get('genome_key'))
        self.fasta = self.genome_data.get('fasta')
        os.makedirs(os.path.join(kwargs.get('output'), "bed_tmp"), exist_ok=True)
        self.runsheet_data = make_scrunsheet(folder=kwargs.get('folder'), output=kwargs.get('output'), \
            sample_flag = kwargs.get('sample_flag'), genome_key = kwargs.get('genome_key'), \
            organized_by="file", strsplit= kwargs.get('strsplit'), ext=kwargs.get('ext'),  r1_char=kwargs.get('r1_char'), \
            r2_char=kwargs.get('r2_char'), fname=kwargs.get('fname'), no_pipe=kwargs.get('no_pipe'), select = kwargs.get('select'))
        self.environs.array_count = len(self.runsheet_data)
        self.commands = self.scalign_executable()
        self.environs.array_job = True
        self.bash_scripts = self.environs.generate_job(self.commands, self.job)
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

    def scalign_executable(self):
        commandline=""
        command = []
        JOBSTRING = self.id_generator(size=10)
        norm_bowtie_flags='--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant -q --phred33 -I 10 -X 700'
        sam2bed_string = """| samTobed - -o %s %s""" % ('${bed_out}', self.filter_string)
        if self.pipe:
            #commandline = commandline + "sed -i 's/\r//g' array_data.tsv\n"
            commandline = commandline + """export fq1=$(head -n $SLURM_ARRAY_TASK_ID array_data.tsv | tail -n 1 | cut -d ',' -f 1 -)\n"""
            commandline = commandline + """export fq2=$(head -n $SLURM_ARRAY_TASK_ID array_data.tsv | tail -n 1 | cut -d ',' -f 2 -)\n"""
            commandline = commandline + """export bed_out=$(head -n $SLURM_ARRAY_TASK_ID array_data.tsv | tail -n 1 | cut -d ',' -f 3 -)\n"""
            commandline = commandline + """bowtie2 %s -p %s -1 ${fq1} -2 ${fq2} -x %s %s\n""" % (self.bowtie_flags, self.threads, self.fasta, sam2bed_string)
            #commandline = commandline + """rm %s \n""" % ('${bed_out}_tmp')
            #commandline = commandline + """rm err.log errtmp out.log outtmp \n"""
        # else:
        #     commandline = """echo '\n[BOWTIE] Running Bowtie for main alignment... Output:\n'\nbowtie2 %s -p %s -1 %s -2 %s -x %s -S %s\n""" % (self.bowtie_flags, self.threads, fastq1, fastq2, sample['fasta'], sample['sam'])
        #     commandline = commandline + """sleep 20s \nsamtools view -bS %s -o %s\n""" % (sample['sam'], sample['bam']+'US')
        #     commandline = commandline + """sleep 20s \necho '\n[SAMTOOLS]... Sorting bam file %s\n'\nsamtools sort -o %s -O bam -T %s %s\n""" % (sample['bam'], sample['bam'], sample['bam']+'US', sample['bam']+'US')
        #     commandline = commandline + """sleep 20s \necho '\n[SAMTOOLS]... Indexing bam file %s\n'\nsamtools index %s\n""" % (sample['bam'],sample['bam'])
        #     commandline = commandline + """sleep 20s \necho '\n[SAMTOBED] Running samToBed for main alignment... Output:\n'\nsamTobed %s -o %s %s""" % (sample['sam'], sample['bed_out']+'tmp', self.filter_string)
        #     commandline = commandline + """\necho '[SORT] Sorting Bed...\n'\nsort -k1,1 -k2n,2n %s > %s\n""" % (sample['bed_out']+'tmp', sample['bed_out'])
        #     commandline = commandline + """rm %s \n""" % (sample['sam'])
        #     commandline = commandline + """rm %s %s \n""" % (sample['bed_out']+'tmp', sample['bam']+'US')
        command.append(["array_job", commandline])
        return command


class SCMerge(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(SCMerge, self).__init__(*args, **kwargs)
        self.job = "HENIPIPE_SCMERGE"
        self.output = kwargs.get('output')
        self.method = kwargs.get('stringency')
        self.norm = kwargs.get('norm')
        self.fdr_thresh = kwargs.get('fdr_thresh')
        self.environs.array_job = False
        self.environs.dependency_job = kwargs.get('dependency')
        self.runsheet_data = kwargs.get('runsheet_data')
        self.genome_data = load_genomes(GENOMES_JSON).get(kwargs.get('genome_key'))
        [i.update({"bc" : striptobarcode(i.get("bed_out"))}) for i in self.runsheet_data]
        self.commands = self.scmerge_executable()
        self.bash_scripts = self.environs.generate_job(self.commands, self.job)
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

    def scmerge_executable(self):
        command = []
        commandline = ""
        JOBSTRING = self.id_generator(size=10)
        first_ = self.runsheet_data[0]
        others_ = self.runsheet_data[1:]
        outfile = os.path.join(self.output,"fragments.tsv")
        outbedg = os.path.join(self.output,"fragments.bg")
        outbw = os.path.join(self.output,"fragments.bw")
        outseacr = os.path.join(self.output,"fragments_seacr_peaks")
        chromsizes = self.genome_data.get('genome_sizes')
        commandline = commandline + """mv outarr samToBed.log\nmv errarr bowtie2.log\n"""
        commandline = commandline + """cat %s | awk -v var=%s '{print $1"\t"$2"\t"$3"\t"var}' > merged.bed\n""" % (first_.get("bed_out"), first_.get("bc"))
        #commandline = commandline + """export beds=$(array_data.tsv | tail -n 1 | cut -d ',' -f 1 -)\n"""
        for other in others_:
            commandline = commandline + """cat %s | awk -v var=%s '{print $1"\t"$2"\t"$3"\t"var}' >> merged.bed\n""" % (other.get("bed_out"), other.get("bc"))
        commandline = commandline + """sort -k1,1 -k2n,2n merged.bed | uniq -c  | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' > %s\n""" % (outfile)
        commandline = commandline + """rm merged.bed\nrm -R %s\n""" % (os.path.join(self.output, "bed_tmp"))
        commandline = commandline + """bedtools genomecov -bg -i %s -g %s > %s\n""" % (outfile, chromsizes, outbedg)
        commandline = commandline + """bgzip %s\n""" % (outfile)
        commandline = commandline + """tabix -p bed %s\n""" % (outfile+".gz")
        commandline = commandline + """bedGraphToBigWig %s %s %s\n""" % (outbedg, chromsizes, outbw)
        commandline = commandline + """bash %s %s %s %s %s %s\n""" % (SEACR_SCRIPT, outbedg, self.fdr_thresh, self.norm, "stringent", outseacr)
        commandline = commandline + """bash %s %s %s %s %s %s\n""" % (SEACR_SCRIPT, outbedg, self.fdr_thresh, self.norm, "relaxed", outseacr)
        #commandline = commandline + """rm %s\n""" % (outbedg)
        command.append(['Dependency_job', commandline])
        return command

class Align(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(Align, self).__init__(*args, **kwargs)
        self.bowtie_flags = kwargs.get('bowtie_flags')
        self.job = "HENIPIPE_ALIGN"
        self.pipe = not kwargs.get('no_pipe')
        self.filter_string = self.make_filter_string(kwargs.get('filter')[0], kwargs.get('filter')[1])
        self.norm_method = kwargs.get('norm_method')
        self.commands = self.align_executable()
        self.bash_scripts = self.environs.generate_job(self.commands, self.job)
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
            norm_bowtie_flags='--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant -q --phred33 -I 10 -X 700'
            sam2bed_string = """| samTobed - -o %s %s""" % (sample['bed_out']+'tmp', self.filter_string)
            if self.pipe:
                commandline = """echo '\n[BOWTIE] Running Bowtie piped to samTobed.py for main alignment... Output:\n'\nbowtie2 %s -p %s -1 %s -2 %s -x %s %s\n""" % (self.bowtie_flags, self.threads, fastq1, fastq2, sample['fasta'], sam2bed_string)
                commandline = commandline + """\necho 'Sorting Bed...\n'\nsort -k1,1 -k2n,2n %s > %s\n""" % (sample['bed_out']+'tmp', sample['bed_out'])
                commandline = commandline + """rm %s \n""" % (sample['bed_out']+'tmp')
            else:
                commandline = """echo '\n[BOWTIE] Running Bowtie for main alignment... Output:\n'\nbowtie2 %s -p %s -1 %s -2 %s -x %s -S %s\n""" % (self.bowtie_flags, self.threads, fastq1, fastq2, sample['fasta'], sample['sam'])
                commandline = commandline + """sleep 20s \nsamtools view -bS %s -o %s\n""" % (sample['sam'], sample['bam']+'US')
                commandline = commandline + """sleep 20s \necho '\n[SAMTOOLS]... Sorting bam file %s\n'\nsamtools sort -o %s -O bam -T %s %s\n""" % (sample['bam'], sample['bam'], sample['bam']+'US', sample['bam']+'US')
                commandline = commandline + """sleep 20s \necho '\n[SAMTOOLS]... Indexing bam file %s\n'\nsamtools index %s\n""" % (sample['bam'],sample['bam'])
                commandline = commandline + """sleep 20s \necho '\n[SAMTOBED] Running samToBed for main alignment... Output:\n'\nsamTobed %s -o %s %s""" % (sample['sam'], sample['bed_out']+'tmp', self.filter_string)
                commandline = commandline + """\necho '[SORT] Sorting Bed...\n'\nsort -k1,1 -k2n,2n %s > %s\n""" % (sample['bed_out']+'tmp', sample['bed_out'])
                commandline = commandline + """rm %s \n""" % (sample['sam'])
                commandline = commandline + """rm %s %s \n""" % (sample['bed_out']+'tmp', sample['bam']+'US')
            if self.norm_method == "spike_in":
                commandline = commandline + """echo '\n[BOWTIE] Running Bowtie piped to samTobed.py for spikein... Output:\n'\nbowtie2 %s -p 4 -1 %s -2 %s -x %s | samTobed - -o %s""" % (norm_bowtie_flags, fastq1, fastq2, sample['spikein_fasta'], sample['spikein_bed_out']+'tmp')
                commandline = commandline + """\necho '[SORT] Sorting Bed for spikein...\n'\nsort -k1,1 -k2n,2n %s > %s\n""" % (sample['spikein_bed_out']+'tmp', sample['spikein_bed_out'])
                commandline = commandline + """rm %s \n""" % (sample['spikein_bed_out']+'tmp')
            command.append([sample['sample'], commandline])
        return command







class Scale(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(Scale, self).__init__(*args, **kwargs)
        self.job = "HENIPIPE_SCALE"
        norm_method = kwargs.get('norm_method')
        #self.processor_line = self.norm_processor_line()
        self.norm_values = self.get_norm_values(method = norm_method)
        self.commands = self.norm_executable()
        self.bash_scripts = self.environs.generate_job(self.commands, self.job)
    def __call__():
        pass

    def get_norm_values(self, method):
        if method=="read_count":
            for sample in self.runsheet_data:
                count = 0
                # for line in open(sample['bed_out']).xreadlines(): count += 1
                for line in open(sample['bed_out']): count += 1
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
                for line in open(sample['bed_out']): count += 1
                ncount = 0
                for line in open(sample['spikein_bed_out']): ncount += 1
                sample['scale_factor'] = (1/float(ncount))
        if method=="none":
            for sample in self.runsheet_data:
                sample['scale_factor'] = 1
        return

    def norm_executable(self):
        command = []
        for sample in self.runsheet_data:
            JOBSTRING = self.id_generator(size=10)
            commandline = """genomeCoverageBed -bg -i %s -g %s -scale %s -trackline | pyWriter - %s\n""" % (sample['bed_out'], sample['genome_sizes'], sample['scale_factor'], sample['bedgraph'])
            command.append([sample['sample'], commandline])
        return command


class Bigwig(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(Bigwig, self).__init__(*args, **kwargs)
        self.job = "HENIPIPE_BIGWIG"
        self.commands = self.bw_executable()
        self.bash_scripts = self.environs.generate_job(self.commands, self.job)
    def __call__():
        pass

    def bw_executable(self):
        command = []
        for sample in self.runsheet_data:
            JOBSTRING = self.id_generator(size=10)
            commandline = """echo '\n[BIGWIG] Making Bigwig  %s from %s'""" % (sample['bigwig'], sample['bedgraph'])
            commandline = commandline + """\nbedGraphToBigWig %s %s %s\n""" % (sample['bedgraph'], sample['genome_sizes'], sample['bigwig'])
            command.append([sample['sample'], commandline])
        return command




class SEACR(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(SEACR, self).__init__(*args, **kwargs)
        self.job = "HENIPIPE_SEACR"
        self.method = kwargs.get('stringency')
        self.norm = kwargs.get('norm')
        self.runsheet_data = self.SEACR_match()
        self.processor_line = self.SEACR_processor_line()
        self.commands = self.SEACR_executable()
        self.bash_scripts = self.environs.generate_job(self.commands, self.job)
    def __call__():
        pass

    def SEACR_match(self):
        sk = [i.get('SEACR_key') for i in self.runsheet_data]
        controls_b = [bool(re.search(r'._CONTROL$', i)) for i in sk]
        controls = list(compress(self.runsheet_data, controls_b))
        samples_b = [not i for i in controls_b]
        samples = list(compress(self.runsheet_data, samples_b))
        for sample in samples:
            control_name = sample.get('SEACR_key')+"_CONTROL"
            try:
                control_bed = next(item for item in controls if item["SEACR_key"] == control_name).get('bedgraph')
            except StopIteration:
                raise ValueError("Could not find matching files, make sure it is in the runsheet and/or your aren't selecting it out with the select flag")
            sample.update( {'SEACR_in' : sample.get('bedgraph')})
            sample.update( {'SEACR_control' : control_bed})
        return samples

    def SEACR_executable(self):
        commandline=""
        command = []
        for sample in self.runsheet_data:
            JOBSTRING = self.id_generator(size=10)
            commandline = """echo '\n[SEACR] Running SEACR... Output:\n'\nbash %s %s %s %s %s %s\n""" % (SEACR_SCRIPT, sample['SEACR_in'], sample['SEACR_control'], self.norm, self.method, sample['SEACR_out'])
            command.append([sample['sample'], commandline])
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
        self.to_process = self.Merge_match()
        self.processor_line = self.Merge_processor_line()
        self.commands = self.Merge_executable()
        self.bash_scripts = self.environs.generate_job(self.commands, self.job)
    def __call__():
        pass

    def Merge_match(self):
        key_data = [i.get("MERGE_key") for i in self.runsheet_data]
        bg_data = [i.get("bedgraph") for i in self.runsheet_data]
        merge_dict = dict.fromkeys(key_data, "NotFound")
        samples = []
        for key in merge_dict.keys():
            samples.append({    "sample" : key,
                                "files_to_merge": list(compress(bg_data, is_in(key, key_data))),})
        return(samples)


    def Merge_executable(self):
        commandline=""
        command = []
        for item in self.to_process:
            seperator = ' '
            nfiles = len(item.get("files_to_merge"))
            bedgraph_line = seperator.join(item.get("files_to_merge"))
            bedgraph_out=str(os.path.join(self.out, item.get("sample")))+"_merged.bedgraph"
            JOBSTRING = self.id_generator(size=10)
            #commandline = """echo '\n[MERGE] Merging bedgraphs:\n%s'\nbedtools unionbedg -i %s | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $0"\t"sum/(NF-4+1); }' > %s\nsleep 10\ncut -d$'\t' -f1-3,%s %s > %s\nrm %s\n""" % (bedgraph_line, bedgraph_line, bedgraph_out+'temp', (nfiles+4), bedgraph_out+'temp', bedgraph_out, bedgraph_out+'temp')
            commandline = """echo '\n[MERGE] Merging bedgraphs:\n%s'\nbedtools unionbedg -i %s | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $0"\t"sum/(NF-4+1); }' > %s\n""" % (bedgraph_line, bedgraph_line, bedgraph_out)
            command.append([item['sample'], commandline])
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
        #self.merged = kwargs.get('merged')
        self.out = kwargs.get('out')
        self.norm = kwargs.get('norm')
        self.runsheet_data = self.MACS2_match()
        self.processor_line = self.MACS2_processor_line()
        self.commands = self.MACS2_executable()
        self.script = self.environs.generate_job(self,commands, self.job)
    def __call__():
        pass

    def MACS2_match(self):
        desired_samples = self.runsheet_data
        #desired_samples = parsed_runsheet
        sample_key = [i.get("sample") for i in desired_samples]
        biomatch_data = [i.get("MACS2_key") for i in desired_samples]
        abmatch_data = [i.get("SEACR_key") for i in desired_samples]
        unique_keys = unique(sample_key)
        run_list = []
        for key in unique_keys:
            #find out if file is bio sample or control or ab sample or control by searching lists of the two keys
            #print(key)
            biomatch_key = [biomatch_data[i] for i in which(key, sample_key)]
            is_biomatch_control = [bool(re.search(r'._CONTROL$', i)) for i in biomatch_key]
            abmatch_key = [abmatch_data[i] for i in which(key, sample_key)]
            is_abmatch_control = [bool(re.search(r'._CONTROL$', i)) for i in abmatch_key]
            is_biomatch_control = all_the_same(is_biomatch_control)
            is_abmatch_control = all_the_same(is_abmatch_control)
            if type(is_abmatch_control) is str:
                raise ValueError("Some discrepency between merge_key and MACS2_key ")
            if type(is_biomatch_control) is str:
                raise ValueError("Some discrepency between merge_key and MACS2_key ")
            if not is_abmatch_control and not is_biomatch_control:
                #treatment bed is just key
                treatment_bed = get_key_from_dict_list(desired_samples, {"sample":key}, 'bed_out')
                #for all non-controls we will output a list of relevant files to process using macs2
                #first find biomatch control
                biomatch_control_key = biomatch_key[0]+"_CONTROL"
                #get bed of this control
                control_bed = get_key_from_dict_list(desired_samples, {"sample":sample_key[which(biomatch_control_key, biomatch_data)[0]]}, 'bed_out')
                #get bed of antibody control for treatment
                abmatch_control_key = abmatch_key[0]+"_CONTROL"
                treatment_abcontrol_bed = get_key_from_dict_list(desired_samples, {"sample":sample_key[which(abmatch_control_key, abmatch_data)[0]]}, 'bed_out')
                #get bed of antibody control for control
                control_control_key = get_key_from_dict_list(desired_samples, {"sample":sample_key[which(biomatch_control_key, biomatch_data)[0]]}, 'SEACR_key')+"_CONTROL"
                control_abcontrol_bed = get_key_from_dict_list(desired_samples, {"sample":sample_key[which(control_control_key, abmatch_data)[0]]}, 'bed_out')
                run_list.append({   "MACS2DIFF_treatment": key,
                                    "MACS2CP_treat_sample": treatment_bed,
                                    "MACS2CP_treat_control": treatment_abcontrol_bed,
                                    "MACS2DIFF_control": sample_key[which(biomatch_control_key, biomatch_data)[0]],
                                    "MACS2CP_control_sample": control_bed,
                                    "MACS2CP_control_control": control_abcontrol_bed,
                                    "sample": key})
        return(run_list)






    def MACS2_executable(self):
        commandline=""
        command = []
        for item in self.runsheet_data:
            JOBSTRING = self.id_generator(size=10)
            treat_p = os.path.join(self.out, item["MACS2DIFF_treatment"])
            cont_p = os.path.join(self.out, item["MACS2DIFF_control"])
            if self.cluster=="SLURM":
                modules = """\nsource /app/Lmod/lmod/lmod/init/bash\nmodule load MACS2"""
            else:
                modules = """\n"""
            commandline = """echo '\n[MACS2] Running MACS2 callpeak on sample... Output:\n'\nmacs2 callpeak -B -t %s -c %s -f BEDPE -g hs --nomodel --extsize 147 --outdir %s -n %s\n""" % (item["MACS2CP_treat_sample"], item["MACS2CP_treat_control"], self.out, item["MACS2DIFF_treatment"])
            commandline = commandline + """echo '\n[MACS2] Getting depth of sample... Output:\n'\nstr1=$(egrep "fragments after filtering in control" %s_peaks.xls | cut -d ":" -f2)\n""" % (treat_p)
            commandline = commandline + """echo '\n[MACS2] Running MACS2 callpeak on control... Output:\n'\nmacs2 callpeak -B -t %s -c %s -f BEDPE -g hs --nomodel --extsize 147 --outdir %s -n %s\n""" % (item["MACS2CP_control_sample"], item["MACS2CP_control_control"], self.out, item["MACS2DIFF_control"])
            commandline = commandline + """echo '\n[MACS2] Getting depth of sample... Output:\n'\nstr2=$(egrep "fragments after filtering in control" %s_peaks.xls | cut -d ":" -f2)\n""" % (cont_p)
            commandline = commandline + """echo '\n[MACS2] Running MACS2 bdgdiff... Output:\n'\nmacs2 bdgdiff --t1 %s_treat_pileup.bdg --c1 %s_control_lambda.bdg --t2 %s_treat_pileup.bdg --c2 %s_control_lambda.bdg --d1 $str1 --d2 $str2 -g 60 -l 147 --o-prefix %s --outdir %s\n""" % (treat_p, treat_p, cont_p, cont_p, item["MACS2DIFF_treatment"]+"_v_"+item["MACS2DIFF_control"], self.out)
            commandline = modules + commandline
            command.append(commandline)
        return command

    def MACS2_processor_line(self):
        if self.cluster=="PBS":
            return """select=1:mem=8GB:ncpus=2"""
        if self.cluster=="SLURM":
            return ''



class AUC(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(AUC, self).__init__(*args, **kwargs)
        self.job = "HENIPIPE_AUC"
        self.out = kwargs.get('out')
        self.norm = kwargs.get('norm')
        self.pipe = not kwargs.get('no_pipe')
        self.method = kwargs.get('stringency')
        self.runsheet_data = self.AUC_match()
        self.processor_line = self.AUC_processor_line()
        self.commands = self.AUC_executable()
        self.script = self.environs.generate_job(self,commands, self.job)
    def __call__():
        pass

    def AUC_match(self):
        desired_samples = self.runsheet_data
        sample_key = [i.get("sample") for i in desired_samples]
        biomatch_data = [i.get("MACS2_key") for i in desired_samples]
        abmatch_data = [i.get("SEACR_key") for i in desired_samples]
        unique_keys = unique(sample_key)
        run_list = []
        for key in unique_keys:
            #find out if file is bio sample or control or ab sample or control by searching lists of the two keys
            biomatch_key = [biomatch_data[i] for i in which(key, sample_key)]
            is_biomatch_control = [bool(re.search(r'._CONTROL$', i)) for i in biomatch_key]
            abmatch_key = [abmatch_data[i] for i in which(key, sample_key)]
            is_abmatch_control = [bool(re.search(r'._CONTROL$', i)) for i in abmatch_key]
            is_biomatch_control = all_the_same(is_biomatch_control)
            is_abmatch_control = all_the_same(is_abmatch_control)
            if type(is_abmatch_control) is str:
                raise ValueError("Some discrepency between merge_key and MACS2_key ")
            if type(is_biomatch_control) is str:
                raise ValueError("Some discrepency between merge_key and MACS2_key ")
            if not is_abmatch_control and not is_biomatch_control:
                #treatment bed is just key
                treatment_bed = get_key_from_dict_list(desired_samples, {"sample":key}, 'bedgraph')
                #for all non-controls we will output a list of relevant files to process using macs2
                #first find biomatch control
                biomatch_control_key = biomatch_key[0]+"_CONTROL"
                #get bed of this control
                control_bed = get_key_from_dict_list(desired_samples, {"sample":sample_key[which(biomatch_control_key, biomatch_data)[0]]}, 'bedgraph')
                #get bed of antibody control for treatment
                abmatch_control_key = abmatch_key[0]+"_CONTROL"
                treatment_abcontrol_bed = get_key_from_dict_list(desired_samples, {"sample":sample_key[which(abmatch_control_key, abmatch_data)[0]]}, 'bedgraph')
                #get bed of antibody control for control
                control_control_key = get_key_from_dict_list(desired_samples, {"sample":sample_key[which(biomatch_control_key, biomatch_data)[0]]}, 'SEACR_key')+"_CONTROL"
                control_abcontrol_bed = get_key_from_dict_list(desired_samples, {"sample":sample_key[which(control_control_key, abmatch_data)[0]]}, 'bedgraph')
                run_list.append({   "AUC_DIFF_treatment": key,
                                    "AUC_CP_treat_sample": treatment_bed,
                                    "AUC_CP_treat_control": treatment_abcontrol_bed,
                                    "AUC_DIFF_control": sample_key[which(biomatch_control_key, biomatch_data)[0]],
                                    "AUC_CP_control_sample": control_bed,
                                    "AUC_CP_control_control": control_abcontrol_bed,
                                    "sample": key})
        return(run_list)

    def AUC_executable(self):
        commandline=""
        command = []
        for item in self.runsheet_data:
            self.files_rm = []
            JOBSTRING = self.id_generator(size=10)
            treat_comb = os.path.join(self.out, (item["AUC_DIFF_treatment"]+"_"+item["AUC_DIFF_control"]+"_treatment_combined.bedgraph"))
            cont_comb = os.path.join(self.out, (item["AUC_DIFF_treatment"]+"_"+item["AUC_DIFF_control"]+"_cont_combined.bedgraph"))
            seacr_merge_prefix = os.path.join(self.out, (item["AUC_DIFF_treatment"]+"_"+item["AUC_DIFF_control"]+"_SEACR"))
            peakfile = os.path.join(self.out, (item["AUC_DIFF_treatment"]+"_"+item["AUC_DIFF_control"]+"_SEACR.")+self.method+".bed")
            out_file = os.path.join(self.out, (item["AUC_DIFF_treatment"]+"_"+item["AUC_DIFF_control"]+"_AUC.bed"))
            cp_treat_out = os.path.join(self.out, (os.path.basename(item["AUC_CP_treat_sample"])+'.bz'))
            cp_cont_out = os.path.join(self.out, (os.path.basename(item["AUC_CP_control_sample"])+'.bz'))
            self.files_rm.extend((cp_treat_out, cp_cont_out, treat_comb, cont_comb, cp_treat_out+".tbi", cp_cont_out+".tbi"))
            if self.cluster=="SLURM":
                modules = """\nsource /app/Lmod/lmod/lmod/init/bash\nmodule load bedtools\nmodule load R\nmodule load htslib/1.9\n"""
            else:
                modules = """\nmodule load bedtools\nmodule load R\nmodule load htslib/1.9\n"""
            commandline = """echo '\n[AUC] Merging sample bedgraphs for aggregated peak call...'\nbedtools unionbedg -i %s %s | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum; }' > %s\n""" % (item["AUC_CP_treat_sample"], item["AUC_CP_control_sample"], treat_comb)
            commandline = commandline + """echo '\n[AUC] Merging control bedgraphs for aggregated peak call...'\nbedtools unionbedg -i %s %s | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum; }' > %s\n""" % (item["AUC_CP_treat_control"], item["AUC_CP_control_control"], cont_comb)
            commandline = commandline + """echo '\n[AUC] Running SEACR on merged sample data... Output:\n'\nbash %s %s %s %s %s %s\n""" % (SEACR_SCRIPT, treat_comb, cont_comb, self.norm, self.method, seacr_merge_prefix)
            commandline = commandline + """echo '\n[AUC] Making Tabix files... \n'\nbgzip -c %s > %s\n""" % (item["AUC_CP_treat_sample"], cp_treat_out)
            commandline = commandline + """bgzip -c %s > %s\n""" % (item["AUC_CP_control_sample"], cp_cont_out)
            commandline = commandline + """tabix -S 1 -p bed %s\n""" % (cp_treat_out)
            commandline = commandline + """tabix -S 1 -p bed %s\n""" % (cp_cont_out)
            commandline = commandline + """auc -o %s -p %s %s %s\n""" % (out_file, peakfile, cp_treat_out, cp_cont_out)
            if self.pipe: commandline = commandline + """echo '\n[AUC] Removing files'\nrm {0}\n""".format(" ".join(self.files_rm))
            commandline = modules + commandline
            command.append(commandline)
        return command

    def AUC_processor_line(self):
        if self.cluster=="PBS":
            return """select=1:mem=4GB:ncpus=1"""
        if self.cluster=="SLURM":
            return ''



def get_key_from_dict_list(list_of_dicts, dict_in, key_out):
    #this function takes a list of dicts, returns the value associated with key_out for the element in the list_of_dicts which has the key:value combo given in dict_in
    dict_out = next(item for item in list_of_dicts if item[list(dict_in.keys())[0]] == dict_in.get(list(dict_in.keys())[0]))
    return dict_out.get(key_out)

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


def find_fastqs_by_sample(dir, strsplit, r1_char, r2_char, ext, full_name=True):
    fastqs=[]
    fastq1=fq1=[]
    fastq2=fq2=[]
    sample=[]
    sample_dict=[]
    #sort into read1 and read2
    for file in os.listdir(dir):
        if file.endswith(ext):
            fastqs.extend([file])
            if r1_char in file:
                fastq1.extend([file])
            if r2_char in file:
                fastq2.extend([file])
    #get samples
    sample=[i.split(strsplit)[0] for i in fastq1]
    #get uniqe samples
    sample = list(set(sample))
    for samp in sample:
        #R1matches=list(filter(lambda x: samp+r1_char in x, fastq1)) (taken out because of lane portion of string)
        R1matches=list(filter(lambda x: samp in x, fastq1))
        R1matches=[os.path.join(dir, i) for i in R1matches]
        fastq1s="\t".join(R1matches)
        #R2matches=list(filter(lambda x: samp+r2_char in x, fastq2))  (taken out because of lane portion of string)
        R2matches=list(filter(lambda x: samp in x, fastq2)) 
        R2matches=[os.path.join(dir, i) for i in R2matches]
        fastq2s="\t".join(R2matches)
        sample_dict.extend([{"sample":samp, "fastq1":fastq1s, "fastq2":fastq2s}])
    return(sample_dict)


def find_fastq_mate(dir, r1_char, r2_char, ext, sample_flag=None, full_name=True):
    fastqs=[]
    fastq1=[]
    fastq2=[]
    for file in os.listdir(dir):
        if file.endswith(ext):
            fastqs.extend([file])
            if r1_char in file:
                fastq1.extend([file])
            if r2_char in file:
                fastq2.extend([file])
    fastq1_mate=[]
    for fastq in fastq1:
        #check if present
        put_R2=re.sub(r1_char, r2_char, fastq)
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


def make_scrunsheet(folder, output, sample_flag, genome_key, organized_by, strsplit, ext, r1_char, r2_char, fname, fasta=None, spikein_fasta=None, genome_sizes=None, no_pipe = True, select = None):
    outfile = "array_data.tsv"
    genome_data = load_genomes(GENOMES_JSON).get(genome_key)
    # if(organized_by=="folder"):    
    #     ddir=[x[0] for x in os.walk(folder)]
    #     dat=list(map(find_fastq_mate, ddir))
    #     #good_dat = [i for i in dat if i.get('has_fastq') is True]
    #     dat=list(map(find_fastq_mate, ddir, [r1_char] * len(ddir), [r2_char] * len(ddir), [ext] * len(ddir)))
    #     good_dat = [i for i in good_dat if re.compile(r'.*'+sample_flag).search(i.get('directory_short'))]
    #     for i in good_dat:
    #         i.update({'sample': i.get('directory_short')})
    #     if(len(good_dat))<1:
    #         raise ValueError("No fastqs found")
    #     else:
    #         print("Found "+str(len(good_dat))+" fastq pairs")
    if(organized_by=="file"):
        good_dat=find_fastqs_by_sample(folder, strsplit=strsplit, r1_char=r1_char, r2_char=r2_char, ext = ext)
        if(len(good_dat))<1:
            raise ValueError("No fastqs found")
        else:
            print("Found "+str(len(good_dat))+" fastq pairs")
    for i in good_dat:
        #i.update({"directory_short":i.get('sample')})
        i.update({'bed_out': os.path.join(output, "bed_tmp", i.get('sample')+".bed")})
            #'spikein_bed_out': os.path.join(output, i.get('directory_short')+"_spikein.bed"), \
            #'bedgraph': os.path.join(output, i.get('directory_short')+".bedgraph"), \
            #'fasta': genome_data.get('fasta'), 'genome_sizes':  genome_data.get('genome_sizes')})
        # if no_pipe:
        #     i.update({'sam': os.path.join(output, i.get('directory_short')+".sam"), \
        #         'bam': os.path.join(output, i.get('directory_short')+".bam")})
        keys = ["fastq1", "fastq2", "bed_out"]
        # if no_pipe:
        #     keys.append("sam")
        #     keys.append("bam")
    #print(select)
    if(select is None):
        select=[x for x in range(0, len(good_dat))]
    #print([good_dat for i in sample_flag])
    final_dat = [good_dat[i] for i in select]
    with open(outfile, 'w') as output_file:
        dict_writer = csv.DictWriter(output_file, fieldnames = keys, extrasaction='ignore', lineterminator='\n')
        #dict_writer.writeheader()
        dict_writer.writerows(final_dat)
    return(final_dat)



def make_runsheet(folder, output, sample_flag, genome_key, organized_by, strsplit, ext, r1_char, r2_char, fname, fasta=None, spikein_fasta=None, genome_sizes=None, no_pipe = True):
    genome_data = load_genomes(GENOMES_JSON).get(genome_key)
    if output is None:
        output = os.path.join(os.getcwd())
    if(organized_by=="folder"):    
        ddir=[x[0] for x in os.walk(folder)]
        #dat=list(map(find_fastq_mate, ddir))
        #good_dat = [i for i in dat if i.get('has_fastq') is True]
        dat=list(map(find_fastq_mate, ddir, [r1_char] * len(ddir), [r2_char] * len(ddir), [ext] * len(ddir)))
        good_dat = [i for i in dat if i.get('has_fastq') is True]
        good_dat = [i for i in good_dat if re.compile(r'.*'+sample_flag).search(i.get('directory_short'))]
        for i in good_dat:
            i.update({'sample': i.get('directory_short')})
        if(len(good_dat))<1:
            raise ValueError("No fastqs found")
        else:
            print("Found "+str(len(good_dat))+" fastq pairs")
    if(organized_by=="file"):
        good_dat=find_fastqs_by_sample(folder, strsplit=strsplit, r1_char=r1_char, r2_char=r2_char, ext = ext)
        if(len(good_dat))<1:
            raise ValueError("No fastqs found")
        else:
            print("Found "+str(len(good_dat))+" fastq pairs")
    for i in good_dat:
        i.update({"directory_short":i.get('sample')})
        i.update({'bed_out': os.path.join(output, i.get('directory_short')+".bed"), \
            'spikein_bed_out': os.path.join(output, i.get('directory_short')+"_spikein.bed"), \
            'bedgraph': os.path.join(output, i.get('directory_short')+".bedgraph"), \
            'bigwig': os.path.join(output, i.get('directory_short')+".bigwig"), \
            'SEACR_key': i.get('directory_short'), \
            'MERGE_key': i.get('directory_short'), \
            'SEACR_out': os.path.join(output, i.get('directory_short')+"_SEACR"), \
            'fasta': genome_data.get('fasta'), 'spikein_fasta': genome_data.get('spikein_fasta'), 'genome_sizes':  genome_data.get('genome_sizes')})
        if no_pipe:
            i.update({'sam': os.path.join(output, i.get('directory_short')+".sam"), \
                'bam': os.path.join(output, i.get('directory_short')+".bam")})
        keys = ["sample", "SEACR_key", "MERGE_key", "fasta", "spikein_fasta", "genome_sizes", "fastq1", "fastq2", "bed_out", "spikein_bed_out", "bedgraph", "bigwig", "SEACR_out"]
        if no_pipe:
            keys.append("sam")
            keys.append("bam")
    with open(os.path.join(output, fname), 'w') as output_file:
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
    required_args = ['sample', 'fastq1', 'fastq2', 'bed_out', 'bedgraph', 'genome_sizes', 'bigwig']
    if args.norm_method == "spike_in":
        required_args.extend(['spikein_fasta', 'spikein_bed_out'])
    if args.job == "SEACR":
        required_args.extend(['SEACR_key', 'SEACR_out'])
    if args.job == "MERGE":
        required_args.extend(['MERGE_key'])
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


def change_dirname(input, substitute):
    intpath = str(input).replace(str(os.path.dirname(input)), substitute)
    outpath = pathlib.Path(intpath)
    return(outpath)

def striptobarcode(input, path=True, ext=".bed"):
    if path:
        out = os.path.basename(input)
    else:
        out = input
    out = re.sub(ext, "", out)
    return(out)
