import argparse
import contextlib
import itertools
import glob
import io
import json
import os
import posixpath
import sys
import tarfile
from functools import reduce
from operator import add

import Bio.SeqIO
import Bio.SeqRecord
import Bio.Seq
import tqdm




def load_mibig_metadata(path):
    with tarfile.open(path, "r|gz") as tar:
        metadata = {}
        for entry in iter(tar.next, None):
            if not entry.name.endswith(".json") or entry.name.startswith("."):
                continue
            bgc = json.load(tar.extractfile(entry))
            if "general_params" in bgc:
                data = bgc["general_params"]
            else:
                data = bgc["cluster"]
            metadata[ data["mibig_accession"] ] = data
    return metadata

def load_mibig_records(path):
    with tarfile.open(path, "r:gz") as tar:
        records = {}
        for entry in iter(tar.next, None):
            if not entry.name.endswith(".gbk") or entry.name.startswith("."):
                continue
            bgc = sum(Bio.SeqIO.parse(io.TextIOWrapper(tar.extractfile(entry)), "genbank"), start=Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("")))
            bgc.id = posixpath.splitext(posixpath.basename(entry.name))[0]
            records[bgc.id] = bgc
    return records


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="The GenBank file obtained after annotation by AntiSMASH.")
parser.add_argument("--list", required=True, help="The list of contigs selected at the deduplication stage.")
parser.add_argument("--mibig", required=True, help="The path to MIBiG data archives.")
parser.add_argument("--version", required=True, choices={"1.3", "2.0"}, help="MIBiG version for which to generate embeddings")
parser.add_argument("--output-prefix", required=True, help="The prefix for the output files.")
parser.add_argument("--retire", help="The list of BGCs to retire")
args = parser.parse_args()

with open(args.list) as f:
    selected = {line.strip() for line in f}

if args.retire is not None:
    with open(args.retire) as f:
        retired = { line.strip().split("\t")[0] for line in f if line.strip() }
else:
    retired = {}

mibig_all_meta_set = set()
for tar in glob.glob(os.path.join(args.mibig, "*.json.tar.gz")):
    mibig_all_meta_set.update(load_mibig_metadata(tar))
mibig_all_meta = sorted(set(mibig_all_meta_set))
print("Loaded", len(mibig_all_meta), "unique MIBiG BGCs")

contigs = []
n_skipped = 0
for record in tqdm.tqdm(Bio.SeqIO.parse(args.input, "genbank"), total=len(selected), desc="Loading contigs"):
    if any(feature.type == "region" for feature in record.features):
        n_skipped += 1
        continue
    record.id = record.name = record.annotations["structured_comment"]["antiSMASH-Data"]["Original ID"]
    contigs.append(record)

print("Skipped", n_skipped, "contigs that still contained biosynthetic regions")
contigs.sort(key=lambda r: r.id)
print("Loaded", len(contigs), "contigs")

if len(contigs) < len(mibig_all_meta):
    print(f"Error: not enough contigs for all BGCs ({len(contigs)} contigs, {len(mibig_all_meta)} bgcs)")


pairing = dict(zip(mibig_all_meta, contigs))
print("Generated static pairs to make sure (contigs, BGC) mapping is consistent across MIBiG versions")

print("Loading MIBiG records for version", args.version)
mibig_gbk = load_mibig_records(os.path.join(args.mibig, f"mibig-{args.version}.gbk.tar.gz"))
print("Loading MIBiG metadata for version", args.version)
mibig_meta = load_mibig_metadata(os.path.join(args.mibig, f"mibig-{args.version}.json.tar.gz"))

with contextlib.ExitStack() as ctx:
    gbk = ctx.enter_context(open(f"{args.output_prefix}.gbk", "w"))
    fna = ctx.enter_context(open(f"{args.output_prefix}.fna", "w"))
    tsv = ctx.enter_context(open(f"{args.output_prefix}.clusters.tsv", "w"))

    print("sequence_id", "bgc_id", "start", "end", "type", file=tsv, sep="\t")
    n_skipped = 0
    n_retired = 0

    for bgc_id, bgc_gbk in tqdm.tqdm(mibig_gbk.items(), total=len(mibig_gbk), desc="Embedding"):

        if bgc_id not in mibig_meta:
            n_skipped += 1
            continue

        if bgc_id in retired:
            n_retired += 1
            continue

        contig = pairing[bgc_id]

        contig_cds = [f for f in contig.features if f.type == "CDS"]
        i = len(contig_cds) // 2
        while contig_cds[i].location.end >= contig_cds[i+1].location.start:
            i += 1

        cds_before = contig_cds[i]
        cds_after = contig_cds[i+1]
        assert cds_before.location.end < cds_after.location.start
        insert_loc = (cds_before.location.end + cds_after.location.start) // 2

        before = contig[:insert_loc]
        after = contig[insert_loc:]
        record = before + bgc_gbk + after

        record.id = record.name = contig.id
        record.annotations["molecule_type"] = "DNA"

        Bio.SeqIO.write(record, gbk, "genbank")
        Bio.SeqIO.write(record, fna, "fasta")

        types = ";".join(sorted(ty for ty in mibig_meta[bgc_id]["biosyn_class"] if ty != "Other"))
        if types == "Nucleoside": # fix bad annotation of BGC0000880 in MIBiG 1.3 JSON
            types = "Unknown"
        print(record.id, bgc_id, insert_loc, insert_loc + len(bgc_gbk), types, file=tsv, sep="\t")

    print("Skipped", n_skipped, "GenBank records that were not in metadata")
    print("Forcefully retired", n_retired, "BGCs")
