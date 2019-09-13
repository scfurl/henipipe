#!/usr/bin/env python
import sys

def run_pyWriter():
    # use stdin if it's full
    if not sys.stdin.isatty():
        input_stream = sys.stdin

    # otherwise, read the given filename
    else:
        try:
            input_filename = sys.argv[1]
        except IndexError:
            message = 'need filename as first argument if stdin is not full'
            raise IndexError(message)
        else:
            input_stream = open(input_filename, 'rU')

    #open outputfile
    try:
        output_filename = sys.argv[2]
    except IndexError:
        message = 'need filename as second argument'
        raise IndexError(message)
    else:
        print("[NORM] Normalization complete\n[GENOMECOVERAGEBED] Making bedgraph...\n")
        output = open(output_filename, 'w')
        i = 0
        for line in input_stream:
            output.writelines(str(line)) # do something useful with each line
            i += 1
        print("\n[GENOMECOVERAGEBED] Wrote %s lines to bedgraph %s...\n" % (i, output_filename))
        output.close()
if __name__ == '__main__':
    run_pyWriter()