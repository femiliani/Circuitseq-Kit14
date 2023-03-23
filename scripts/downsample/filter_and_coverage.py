# This script takes an input file, downsamples it using the length to make it more likely to work with flye

#usage python filter_and_coverage.py --infile X --outfile X --nxx 70 --coverage 60

import glob
import gzip
import argparse

def calculate_n70(seq_lengths, cutoff):

    # Sort the sequence lengths in descending order
    seq_lengths.sort(reverse=True)

    # Calculate the total length of the sequences
    total_length = sum(seq_lengths)

    # Calculate the N value
    n70 = 0
    running_total = 0
    for length in seq_lengths:
        running_total += length
        if running_total >= total_length * cutoff/100:
            n70 = length
            break

    return(n70)

def filter_fq(file, outfile):
    # if file is gzipped
     
    try:
        with gzip.open(file, 'rb') as f:
            lengths = []
            for header, seq, plus, qual in zip(*[f]*4):
                lengths.append(len(seq.rstrip()))
            N70 = calculate_n70(lengths, cutoff)
        
        count = 0
        with gzip.open(file, 'rb') as f, gzip.open(outfile, 'wb') as out:
            for header, seq, plus, qual in zip(*[f]*4):
                if len(seq.rstrip()) > round(N70*.9):

                    out.write(header+seq+plus+qual)
                    count +=1 
                if count > coverage:
                    break
    # if it isnt
    except OSError:

        with open(file, 'r') as f:
            lengths = []
            for header, seq, plus, qual in zip(*[f]*4):
                lengths.append(len(seq.rstrip()))
            N70 =calculate_n70(lengths, cutoff)
        
        count = 0
        with open(file, 'r') as f, open(outfile, 'w') as out:
            for header, seq, plus, qual in zip(*[f]*4):
                if len(seq.rstrip()) > round(N70*.9):

                    out.write(header+seq+plus+qual)
                    count +=1 
                if count > coverage:
                    break
            
    return(print(f'Returned {coverage}x coverage at N{cutoff} = {N70}bp and min read length = {round(N70*.9)}bp'))


parser = argparse.ArgumentParser(description='Take input fastq, output 40x coverage of N70 for Flye')

# Add a required input file argument
parser.add_argument('--infile', metavar='INFILE', type=str,
                    help='input file path')

# Add a required output file argument
parser.add_argument('--outfile', metavar='OUTFILE', type=str,
                    help='output file path')

# Add a required Nx argument
parser.add_argument('--nxx', metavar='NXX', type=float,
                    help='length cut off, eg 50 for N50')

# Add a required output file argument
parser.add_argument('--coverage', metavar='COVERAGE', type=float,
                    help='desired coverage, eg 40 for 40x coverage ')

# Parse the command line arguments
args = parser.parse_args()



# Save the input and output file paths as variables
infile = args.infile
outfile = args.outfile
cutoff = args.nxx
coverage = args.coverage


filter_fq(infile, outfile)
