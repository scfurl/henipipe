#!/usr/bin/env python

#AUC finder
#this script will find AUC measurements using ranges taken from a peak_file and interrogating a number of target files
#requires tabix

#whole_peak and narrow_peak
import argparse
import os
from subprocess import Popen, PIPE
import time




class AUC:
    def __init__(self, *args, **kwargs):
        self.peak_file = kwargs.get('peak_file')
        self.targets = kwargs.get('targets')
        self.file_out = kwargs.get('file_out')
        self.header = kwargs.get('header')
        self.peak_count=0

    def __call__():
        pass

    def address_ok(self, address):
        first_split = address.split(":")
        locs = first_split[1].split("-")
        if any([not s for s in first_split]) or any([not s for s in locs]):
            return False
        else:
            return True

    def peak_line(self, peak):
        address = peak.split("\t")[:3]
        wide_peak = "{0:s}:{1:s}-{2:s}".format(*address)
        narrow_peak = peak.split("\t")[5].split()[0]
        #i = targets[0]
        #print(wide_peak)
        #print(narrow_peak)
        if self.address_ok(wide_peak) and self.address_ok(narrow_peak):
            wide_values=[]
            narrow_values=[]
            for i in self.targets:
                proc = Popen("tabix {0:s} {1:s} | awk '{{s+=$4}} END {{print s}}'".format(i, wide_peak), shell = True, stdin = PIPE, stdout=PIPE, stderr = PIPE)
                out, err = proc.communicate()
                try:
                    wide_values.append(out.decode('UTF-8').split()[0])
                except IndexError:
                    wide_values.append("NF")
                proc = Popen("tabix {0:s} {1:s} | awk '{{s+=$4}} END {{print s}}'".format(i, narrow_peak), shell = True, stdin = PIPE, stdout=PIPE, stderr = PIPE)
                out, err = proc.communicate()
                try:
                    narrow_values.append(out.decode('UTF-8').split()[0])
                except IndexError:
                    narrow_values.append("NF")
            return "\t".join(["\t".join(address), ("\t".join(wide_values)), ("\t".join(narrow_values)), narrow_peak])+"\n"
        if self.address_ok(wide_peak) and not self.address_ok(narrow_peak):
            wide_values=[]
            print("Narrow peak error found on line {0} of peak file".format(self.peak_count))
            for i in self.targets:
                proc = Popen("tabix {0:s} {1:s} | awk '{{s+=$4}} END {{print s}}'".format(i, wide_peak), shell = True, stdin = PIPE, stdout=PIPE, stderr = PIPE)
                out, err = proc.communicate()
                try:
                    wide_values.append(out.decode('UTF-8').split()[0])
                except IndexError:
                    wide_values.append("NF")
            return "\t".join(["\t".join(address), ("\t".join(wide_values)), ("\t".join(["NF"]*len(self.targets))), narrow_peak])+"\n"
        if not self.address_ok(wide_peak):
            print("Wide peak error found on line {0} of peak file".format(self.peak_count))
            return None

    def calculate_AUCs(self, start=1, end=1e6):
        """
        peak_file = open(args.peak_file, 'r')
        targets = args.targets
        file_out = open(args.output, 'w')
        """
        start_time = time.time()
        peak_file = open(self.peak_file, 'r')
        file_out = open(self.file_out, 'w')
        if self.header:
            target_bn_wide = "\t".join([os.path.basename(i)+"_wide" for i in self.targets])
            target_bn_narrow = "\t".join([os.path.basename(i)+"_narrow" for i in self.targets])
            header = "seq\tbegin\tend\t"+target_bn_wide+"\t"+target_bn_narrow+"\tnarrow_coord\n"
            file_out.writelines(header)
        iterator = iter(peak_file)
        done_looping = False
        while not done_looping:
            try:
                peak = next(iterator)
                self.peak_count += 1
                #print(self.peak_count)
                if self.peak_count >= end:
                    done_looping = True #2265 error
                    print("\n[AUC] Output: \n Processed {0} peaks\n".format(self.peak_count))
            except StopIteration:
                peak_file.close()
                file_out.close()
                print("\n[AUC] Output: \n Processed {0} peaks".format(self.peak_count))
                done_looping = True
                end_time = time.time()
                print("\n[AUC] Output: \n Execution took {0} seconds\n".format(end_time - start_time))
            else:
                #line_out = peak_line(peak, args.targets)
                if self.peak_count >= start:
                    line_out = self.peak_line(peak)
                    if line_out is not None: file_out.writelines(line_out)
        return


def run_auc():
    parser = argparse.ArgumentParser('A script that will find AUC measurements using ranges taken from a peak_file and interrogating a number of target files')
    parser.add_argument('--peak_file', '-p', type=str, help='A peak file (only SEACR supported)')
    parser.add_argument('--output', '-o', type=str, default=".", help='Write files to this location')
    parser.add_argument('--keepfiles', '-k', action ='store_true', default=False,  help='Write files to this location')
    parser.add_argument('--noheader', '-nh', action ='store_true', default=False,  help='Do not include a header')
    parser.add_argument('targets', nargs='*')
    args = parser.parse_args()
    """
    call = 'auc -o fc/4N1_H3K27me3_4G1_H3K27me3_AUC.bed -p /active/furlan_s/Data/CNR/190801_CNRNotch/henipipe_150/fc/4N1_H3K27me3_4G1_H3K27me3_SEACR.stringent.bed /active/furlan_s/Data/CNR/190801_CNRNotch/henipipe_150/fc/4N1_H3K27me3.bedgraph.bz /active/furlan_s/Data/CNR/190801_CNRNotch/henipipe_150/fc/4G1_H3K27me3.bedgraph.bz'
    args = parser.parse_args(call.split(" ")[1:])
    """

    if os.path.isabs(args.output) is False:
        args.output = os.path.abspath(args.output)

    AUCjob = AUC(peak_file = args.peak_file, targets = args.targets, file_out = args.output, header = not args.noheader)
    # AUCjob.targets
    # AUCjob.peak_file
    AUCjob = AUCjob.calculate_AUCs()

if __name__ == '__main__':
    run_auc()




