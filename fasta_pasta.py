# -*- coding: utf-8 -*-
"""
Extract fasta sequences from a multifasta file using list of sequence ids

@Authors: Michael Vacher
@Date: 2015.05.05
"""

import argparse
import sys
import csv
import time
import re


def get_args(argv):
    desc = '''Extract fasta sequences from a multifasta file using list of sequence ids'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-o','--output_file', help='output file (fasta)', required=False)
    parser.add_argument('-f','--fasta', help='Main fasta file', required=True)
    parser.add_argument('-i','--id_file', help='Sequence id (if multiple columns, the id must be in the first one)', required=True)
    args = parser.parse_args()
    return args

def get_seq_ids(args):
    '''
    Read the id file and return a dictionary of id
    values are boolean indicating if the seq has been processed
    :param args: mains args
    :return: list of id (string)
    '''
    seq_ids = {}
    with open(args.id_file) as tsv:
        for line in csv.reader(tsv, delimiter="\t"):
            id = line[0]
            #skip comments
            if id.startswith("#"):
                continue
            if seq_ids.has_key(id):
                continue
                #print("Id %s is duplicated in %s"%(id, args.id_file))
            seq_ids[id] = False
    return seq_ids

def read_fasta(fp):
    '''
    :param fp: File handler
    :return:
    '''
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def write_fasta_record(fp, title, data, wrap=60):
    '''
    :param fp: file handler (write)
    :param title: sequence id (without >)
    :param data: sequence data (no space no new line
    :param wrap: length of the lines
    :return:
    '''
    """Write a single Fasta record to the file."""
    assert "\n" not in title
    assert "\r" not in title
    fp.write(">%s\n" % title)


    assert "\n" not in data
    assert "\r" not in data
    if wrap:
        for i in range(0, len(data), wrap):
            fp.write(data[i:i + wrap] + "\n")
    else:
        fp.write(data + "\n")

def extract_sequences(args, seq_ids):
    '''
    :param args: mains args
    :param seq_ids: dictionary of sequence id
    :return:
    '''
    main_fasta_ids = {}
    with open(args.output_file, 'w') as output_file:
        with open(args.fasta) as fp:
            for full_id, seq in read_fasta(fp):
                id =  re.sub('^>', '', full_id)
                if main_fasta_ids.has_key(id):
                    print("Id %s is duplicated in %s"%(id, args.fasta))

                main_fasta_ids[id] = False
                if id in seq_ids.keys():
                    write_fasta_record(output_file, id, seq)
                    seq_ids[id] = True
                    main_fasta_ids[id] = True

    return main_fasta_ids, seq_ids

if __name__ == "__main__":
    # Get args
    args = get_args(sys.argv[1:])

    start_time = time.time()
    seq_ids = get_seq_ids(args)
    main_ids, seq_ids = extract_sequences(args, seq_ids)

   # print some final stats
    print "\n","#"*60,"\n"
    print "%s: \tTotal number of sequences %i\t(%i processed)"%(args.fasta, len(main_ids.keys()), sum(main_ids.values()))
    print "%s: \tTotal number of sequences %i\t(%i processed)"%(args.id_file, len(seq_ids.keys()), sum(seq_ids.values()))
    print "\t\t-> %i sequences saved to %s"%(sum(main_ids.values())+sum(seq_ids.values()),args.id_file)
    print 'Total processing time: %.3fs'%(time.time()-start_time)
    print "\n","#"*60,"\n"