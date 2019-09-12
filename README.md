
#henipipe

A python wrapper for processing of sequencing data using CutnRun or CutnTag

##Installation

```bash
pip install

```

##Usage

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

##Examples

Say your fastqs live here
```bash
cd /data/fastq
```

Then
```bash
cd ..
mkdir henipipe
cd henipipe
henipipe MAKERUNSHEET -fq ../fastq -sf MySamplesStartWithThisString -o henipipe
henipipe ALIGN -r runsheet.csv
henipipe NORM -r runsheet.csv
henipipe SEACR -r runsheet.csv
```