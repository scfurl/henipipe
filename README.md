# henipipe
==========  

A python wrapper for processing of sequencing data generated using CutnRun or CutnTag (developed by the Henikoff lab FHCRC)

## Requirements

1. Python > 3.5
2. Computing cluster with PBS or SLURM
3. Modules installed for python, bowtie2, samtools, bedtools, R, SEACR

## Installation

Load python 3 then install henipipe
```bash
module load python
pip install henipipe --user
```

If you run into errors that look like this...
```bash
Retrying (Retry(total=4, connect=None, read=None, redirect=None, status=None)) after connection broken\n by 'SSLError(SSLError(1, u'[SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed (_ssl.c:618)\n'),)': /simple/henipipe/
```

Try rerunning with trusted host options
```bash
pip install henipipe --trusted-host pypi.org --trusted-host files.pythonhosted.org --user
```

Unfortunately, current version of henipipe requires PATH prepending.  So finally you will need to run something like this below to run henipipe directly from the command line.
```bash
python3 -m site &> /dev/null && PATH="$PATH:`python3 -m site --user-base`/bin"
```

You should then be able to test installation by calling henipipe.  You should see the help screen displayed.

```bash
henipipe
```

## Usage

```bash
henipipe [-h] [--sample_flag SAMPLE_FLAG]
                                      [--fastq_folder FASTQ_FOLDER]
                                      [--filter_high FILTER_HIGH]
                                      [--filter_low FILTER_LOW]
                                      [--output OUTPUT] [--runsheet RUNSHEET]
                                      [--log_prefix LOG_PREFIX]
                                      [--select SELECT] [--debug]
                                      [--bowtie_flags BOWTIE_FLAGS]
                                      [--cluster {PBS,SLURM}]
                                      [--norm_method {coverage,read_count,spike_in}]
                                      [--user USER] [--SEACR_norm {non,norm}]
                                      [--SEACR_stringency {stringent,relaxed}]
                                      [--verbose]
                                      {MAKERUNSHEET,ALIGN,NORM,SEACR}

positional arguments:
  {MAKERUNSHEET,ALIGN,NORM,SEACR}
                        a required string denoting segment of pipeline to run.
                        1) "MAKERUNSHEET" - to parse a folder of fastqs; 2)
                        "ALIGN" - to perform alignment using bowtie and output
                        bed files; 3) "NORM" - to normalize data to reference
                        (spike in); 4) "SEACR" - to perform SEACR.

optional arguments:
  -h, --help            show this help message and exit
  --sample_flag SAMPLE_FLAG, -sf SAMPLE_FLAG
                        FOR MAKERUNSHEET only string to identify samples of
                        interest in a fastq folder
  --fastq_folder FASTQ_FOLDER, -fq FASTQ_FOLDER
                        For MAKERUNSHEET only: Pathname of fastq folder (files
                        must be organized in folders named by sample)
  --filter_high FILTER_HIGH, -fh FILTER_HIGH
                        For ALIGN only: upper limit of fragment size to
                        exclude, defaults is no upper limit. OPTIONAL
  --filter_low FILTER_LOW, -fl FILTER_LOW
                        For ALIGN only: lower limit of fragment size to
                        exclude, defaults is no lower limit. OPTIONAL
  --output OUTPUT, -o OUTPUT
                        For MAKERUNSHEET only: Pathname to write runsheet.csv
                        file (folder must exist already!!), Defaults to
                        current directory
  --runsheet RUNSHEET, -r RUNSHEET
                        tab-delim file with sample fields as defined in the
                        script. - REQUIRED for all jobs except MAKERUNSHEET
  --log_prefix LOG_PREFIX, -l LOG_PREFIX
                        Prefix specifying log files for henipipe output from
                        henipipe calls. OPTIONAL
  --select SELECT, -s SELECT
                        To only run the selected row in the runsheet, OPTIONAL
  --debug, -d           To print commands (For testing flow). OPTIONAL
  --bowtie_flags BOWTIE_FLAGS, -b BOWTIE_FLAGS
                        For ALIGN: bowtie flags, OPTIONAL
  --cluster {PBS,SLURM}, -c {PBS,SLURM}
                        Cluster software. OPTIONAL Currently supported: PBS
                        and SLURM
  --norm_method {coverage,read_count,spike_in}, -n {coverage,read_count,spike_in}
                        For ALIGN and NORM: Normalization method, by
                        "read_count", "coverage", or "spike_in". If method is
                        "spike_in", HeniPipe will align to the spike_in
                        reference genome provided in runsheet. OPTIONAL
  --user USER, -u USER  user for submitting jobs - defaults to username.
                        OPTIONAL
  --SEACR_norm {non,norm}, -Sn {non,norm}
                        For SEACR: Normalization method; default is
                        "non"-normalized, select "norm" to normalize using
                        SEACR. OPTIONAL
  --SEACR_stringency {stringent,relaxed}, -Ss {stringent,relaxed}
                        FOR SEACR: Default will run as "stringent", other
                        option is "relaxed". OPTIONAL
  --verbose, -v         Run with some additional ouput - not much though...
                        OPTIONAL
```


## Runsheet

The runsheet is the brains of the henipipe workflow.  You can make a runsheet using the MAKERUNSHEET command.  This command will parse a directory of fastq folder (specified using the -fq flag; fastq files should be organized in subfolders named by sample) and will find fastq mates (R1 and R2).  There are optional for selecting only folders that contain a specific string (using the -sf flag).  *After generation of a runsheet (csv file), you should take a look at it in Excel or Numbers to make sure things look okay...*  Here are the columns that you can include.  Order is irrelevant.  Column names (headers) are suggested.

### Example Runsheet 


| sample | fasta | spikein_fasta | fastq1 | fastq2 | bed_out | spikein_bed_out | genome_sizes | bedgraph |  SEACR_key  | SEACR_out |
|--------|-------|---------------|--------|--------|---------|-----------------|--------------|----------|-------------|-----------|
|  mys1  |  path |      path     |  path  |  path  |   path  |       path      |     path     |   path   |     4JS     |   path    |
|  mys2  |  path |      path     |  path  |  path  |   path  |       path      |     path     |   path   | 4JS_CONTROL |   path    |


* 'sample' name of the sample REQUIRED.  
* 'fasta' location of the Bowtie2 indexed fasta file REQUIRED.  
* 'spikein_fasta' location of the Bowtie2 indexed fasta file for spike_in normalization OPTIONAL.  
* 'fastq1' a tab seperated string of filenames denoting location of all R1 files for a sample REQUIRED.  
* 'fastq2' a tab seperated string of filenames denoting location of all R2 files for a sample REQUIRED.  
* 'bed_out' name of the location for the aligned and sorted bam file REQUIRED.  
* 'spikein_bed_out' name of the location for the aligned and sorted bam file OPTIONAL.  
* 'genome_sizes' REQUIRED.  
* 'bedgraph' file name of normalized bedgraph REQUIRED.  
* 'SEACR_key' sample key corresponding to sample groups to be run against an IgG (or other) contol.  all samples to be run against a control are given the same name and the control is labeleled with the an additional string underscore + 'CONTROL' (i.e. 4JS_CONTROL) OPTIONAL.  
* 'SEACR_out' file name of SEACR output OPTIONAL.  



## Examples

Say your fastqs live within within subfolders of a folder 'fastq' in the folder 'data'.  So if you were to...
```bash
cd /data/fastq
ls
```
... you'd get a bunch of folders, each of which would be filled with fastqs.  Each folder name should correspond to a sample name.


To run henipipe, do the following...
1. Make a new output directory 'henipipe'.
2. Go into that directory and make a runsheet pointing to the fastq folder ( folder level above.  (henipipe is cool with realtive and absolute pathnames here; but as stated earlier, absolute pathnames are best for the runsheet.)  
3.  Optionally you can only select directories of fastq files that contain the string denoted using the -sf flag.  
4. After inspecting and completing the runsheet, run ALIGN, NORM, and SEACR.  
5. Sit back have a cocktail.

```bash
cd ..
mkdir henipipe
cd henipipe
henipipe MAKERUNSHEET -fq ../fastq -sf MySamplesStartWithThisString -o henipipe
henipipe ALIGN -r runsheet.csv
henipipe NORM -r runsheet.csv
henipipe SEACR -r runsheet.csv
```


##Acknowledgements

Written by Scott Furlan with code inspiration from Andrew Hill's cellwrapper; Henipipe includes a python script sam2bed2.py which takes code from a fantastic sam reader "simplesam" - https://github.com/mdshw5/simplesam