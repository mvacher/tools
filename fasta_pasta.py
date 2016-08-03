# -*- coding: utf-8 -*-
from __future__ import division, print_function
import argparse
import sys
import csv
import time
import re
import os

"""
Extract fasta sequences from a multifasta file :
    - using list of sequence ids
    - using the species id

@Authors: Michael Vacher
@Updated: 2016.08.08
"""


def get_args(argv):
    desc = '''Extract fasta sequences from a multifasta file.
Action:
    'extract': extract fasta sequence using list of sequence ids
    'separate-species': split the fasta file into multiple file (1file/species)
'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-o', '--output_file', help='output file (fasta)',
                        required=False)
    parser.add_argument('--output_directory', help='output directory \
                        (only required with --action=separate-species)',
                        required=False)
    parser.add_argument('-f', '--fasta', help='Main fasta file', required=True)
    parser.add_argument('--action', choices=['extract', 'separate-species'],
                        required=True)
    parser.add_argument('-i', '--id_file', help='Sequence id \
                        (if multiple columns, ids must be in the first one)\
                        Only required with --action=extract',
                        required=False)
    args = parser.parse_args()

    # Check parameters
    if args.action == 'separate-species':
        if not args.output_directory:
            parser.error('--output_directory must be provided when \
--action={}'.format(args.action))
        if args.output_file:
            print("--output_file will be ignored when \
--action={}".format(args.action))
        if not os.path.exists(args.output_directory):
            print("Output directory {} does not exist and will be created"
                  .format(args.output_directory))
            os.makedirs(args.output_directory)
        if args.id_file:
            print("--id_file will be ignored when \
--action={}".format(args.action))
    if args.action == 'extract':
        if args.output_directory:
            print("--output_directory will be ignored when \
--action={}".format(args.action))
        if not args.output_file:
            parser.error('--output_file must be provided when \
--action={}'.format(args.action))
        if not args.id_file:
            parser.error('--id_file must be provided when \
--action={}'.format(args.action))

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
            if id.startswith("#"):  # skip comments
                continue
            if id in seq_ids:
                sys.exit("ERROR: Id {} is duplicated in {}"
                         .format(id, args.id_file))
            seq_ids[id] = False
    return seq_ids


def get_seq_ids_from_fasta(args):
    '''
    Extract ids from a fasta file
    :return: list of id (string)
    '''
    seq_ids = {}
    with open(args.fasta) as fp:
        for full_id, seq in read_fasta(fp):
            id = re.sub('^>', '', full_id)
            if id in seq_ids:
                sys.exit("ERROR: Id {} is duplicated in {}"
                         .format(id, args.id_file))
            seq_ids[id] = False
    return seq_ids


def get_species_ids(args, sequence_ids, delimiter='-'):
    '''
    Convert a list of sequence ids into species ids
    :param sequence_ids: dictionary of sequence ids
    :param delimiter: separator in the sequence id
    :return: dictionary of species id (key: species id, value: list of seq id)
    '''
    species_ids = {}
    for id in seq_ids:
        s = id.split("-")
        assert(len(s) == 2)
        # s[0] = species id; s[1] = secondary sequence id
        if s[0] not in species_ids:
            species_ids[s[0]] = [s[1]]
        else:
            species_ids[s[0]].append(s[1])
    return species_ids


def read_fasta(fp):
    '''
    :param fp: File handler
    :return:
    '''
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))


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
    :param seq_ids: dictionary of sequence id to extract
    :return:
    '''
    main_fasta_ids = {}
    with open(args.output_file, 'w') as output_file:
        with open(args.fasta) as fp:
            for full_id, seq in read_fasta(fp):
                id = re.sub('^>', '', full_id)
                if id in main_fasta_ids:
                    print("Id {} is duplicated in {}".format(id, args.fasta))

                main_fasta_ids[id] = False
                if id in seq_ids.keys():
                    write_fasta_record(output_file, id, seq)
                    seq_ids[id] = True
                    main_fasta_ids[id] = True

    return main_fasta_ids, seq_ids


def separate_sequences_by_species(args, species_ids):
    '''
    Loop over a dictionary of species ids and export sequences from
    each species into separate files
    '''
    for sId, list_of_sec_id in species_ids.items():
        saved_seq = 0
        oFile = os.path.join(args.output_directory, sId + ".fasta")
        print("\tProcessing species {} [{} sequences] - output file: {}"
              .format(sId, len(list_of_sec_id), oFile))

        with open(oFile, 'w') as output_file:
            with open(args.fasta) as fp:
                for full_id, seq in read_fasta(fp):
                    id = re.sub('^>', '', full_id)
                    main_spe_id, sec_id = id.split("-")

                    if main_spe_id == sId:
                        # Check the secondary id is in the list
                        assert sec_id in list_of_sec_id
                        write_fasta_record(output_file, id, seq)
                        saved_seq += 1
        # Verify that the number of saved seq match the initial number
        assert saved_seq == len(list_of_sec_id)


if __name__ == "__main__":
    # Get args
    args = get_args(sys.argv[1:])
    start_time = time.time()

    print("\n", "#"*60, "\n")
    if args.action == 'extract':
        seq_ids = get_seq_ids(args)
        main_ids, seq_ids = extract_sequences(args, seq_ids)
        # print some final stats
        print("{}: \tTotal number of sequences {}\t({} processed)"
              .format(args.fasta, len(main_ids.keys()),
                      sum(main_ids.values())))
        print("{}: \tTotal number of sequences {}\t({} processed)"
              .format(args.id_file, len(seq_ids.keys()),
                      sum(seq_ids.values())))
        print("\t\t-> {} sequences saved to {}"
              .format(sum(main_ids.values())+sum(seq_ids.values()),
                      args.id_file))

    if args.action == 'separate-species':
        print("- Extracting sequence ids from {}".format(args.id_file))
        seq_ids = get_seq_ids_from_fasta(args)
        species_ids = get_species_ids(args, seq_ids)
        separate_sequences_by_species(args, species_ids)
        print("\nTotal number of sequences to extract: {}"
              .format(len(seq_ids.keys())))
        print("Total number of species: {}"
              .format(len(species_ids.keys())))

    print('Total processing time: {}s'.format(time.time()-start_time))
    print("\n", "#"*60, "\n")
