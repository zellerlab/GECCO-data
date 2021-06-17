import argparse
import itertools
import sys
from functools import reduce
from operator import add

import Bio.SeqIO
import Bio.SeqRecord
import tqdm



parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="The GenBank file obtained after annotation by AntiSMASH.")
parser.add_argument("--list", required=True, help="The list of contigs selected at the deduplication stage.")
parser.add_argument("--output", required=True, help="The file where to write the contigs with biosynthetic regions removed.")
args = parser.parse_args()

with open(args.list) as f:
    selected = {line.strip() for line in f}

n_done = 0


with open(args.output, "w") as f:

    for record in tqdm.tqdm(Bio.SeqIO.parse(args.input, "genbank"), total=len(selected)):
        # check the record was processed by AntiSMASH 5.2
        if not record.annotations['structured_comment']['antiSMASH-Data']['Version'].startswith("5.2"):
            raise ValueError(f"record {record.id} was not processed by antiSMASH 5.2")

        # record indices where to clip the sequence
        clip_indices = [0]
        for feature in record.features:
            if feature.type == "region":
                clip_indices.append(feature.location.start)
                clip_indices.append(feature.location.end+1)
        clip_indices.append(len(record))

        # extract subsequences between regions to remove
        subseqs = []
        for i in range(0, len(clip_indices), 2):
            start, end = clip_indices[i:i+2]
            subseqs.append(record.seq[start:end])

        # build new record without regions
        new_record = Bio.SeqRecord.SeqRecord(reduce(add, subseqs))
        new_record.id = record.annotations['structured_comment']['antiSMASH-Data']['Original ID']
        new_record.name = record.annotations['structured_comment']['antiSMASH-Data']['Original ID']
        new_record.annotations['molecule_type'] = 'DNA'

        # write the new record
        Bio.SeqIO.write(new_record, f, "fasta")

        #
        n_done += new_record.id in selected


if n_done < len(selected):
    print("Error: some contigs that should have been selected were not found in the input file")
    sys.exit(1)
