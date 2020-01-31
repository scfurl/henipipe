#!/usr/bin/python
#v0.2

import sys
import argparse


# inFileA = args.sam1
# inFileB = args.sam2

def samComp_writer(inFileA, inFileB, outfile, mode):
    reads = {}
    for line in inFileA:
        #deal with header
        if line[0] != "@": 
            data = line.strip().split('\t')
            name = data[0]
            flag = int(data[1])
            if not flag & 4:
                reads[name] = True

    if mode=="intersect":
        for line in inFileB:
            #deal with header
            if line[0] != "@": 
                data = line.strip().split('\t')
                name = data[0]
                flag = int(data[1])
                if not flag & 4:
                    if name in reads:
                        outfile.writelines('\t'.join(data)+'\n')

    if mode=="unique":
        for line in inFileB:
            #deal with header
            if line[0] != "@":
                data = line.strip().split('\t')
                name = data[0]
                flag = int(data[1])
                if not flag & 4:
                    if name not in reads:
                        outfile.writelines('\t'.join(data)+'\n')
    return


def samComp():
    parser = argparse.ArgumentParser('A script finding either the "intersect" or "unique" reads in two SAM files (for "unique": i.e reads in sam2 not in sam1')
    parser.add_argument('sam1', type=str, default = sys.stdin, help='reference sam')
    parser.add_argument('sam2', type=str, default = sys.stdin, help='reference sam')
    parser.add_argument('--out_sam', '-o', default = "./out.sam", help='SAM output', required=False)
    parser.add_argument('--mode', '-m', type = str, default = "unique", choices=['unique', 'intersect'], help='mode to compare sam files', required=False)
    args = parser.parse_args()
    """
    call='samComp /Users/sfurla/Desktop/smaller.sam /Users/sfurla/Desktop/smaller.sam -m intersect'
    args = parser.parse_args(call.split(" ")[1:])
    """

    #Deal with pipes and open files
    if args.sam1=="-" and args.sam2=="-":
        raise ValueError("Two pipes not supported")
    if args.sam1=="-":
        args.sam1 = sys.stdin
    else:
        args.sam1 = open(args.sam1, 'r')

    if args.sam2=="-":
        args.sam2 = sys.stdin
    else:
        args.sam2 = open(args.sam2, 'r')

    #open write file
    args.out_sam = open(args.out_sam, 'w')

    samComp_writer(args.sam1, args.sam2, args.out_sam, args.mode)

if __name__ == '__main__':
    samComp()




