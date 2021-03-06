import argparse
import glob
import multiprocessing.pool
import os
from functools import partial

import pyfastani
import disjoint_set

import Bio.SeqIO
import tqdm



parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="The input file containing all the contigs to group.")
parser.add_argument("--output", help="The file where to write the selected contigs.")
parser.add_argument("--list", required=True, help="The file where to write the final list of selected contigs.")
parser.add_argument("--6genomes", dest="six", required=True, help="The path to the 6 genomes")
parser.add_argument("--9genomes", dest="nine", required=True, help="The path to the 9 genomes")
parser.add_argument("--cutoff", type=float, default=80, help="The ANI cutoff above which to merge sequences")
parser.add_argument("--total", type=int, default=5000, help="The number of records in the input file")
args = parser.parse_args()

if len(glob.glob(os.path.join(args.six, "*.gbk"))) != 6:
    print("Error: missing sequences in 6 genomes folder")
    exit(1)
if len(glob.glob(os.path.join(args.nine, "*.gbk"))) != 13:
    print("Error: missing sequences in 9 genomes folder")
    exit(1)

sizes = {}
mapper = pyfastani.Mapper()

print("Sketching references...")
for record in tqdm.tqdm(Bio.SeqIO.parse(args.input, "fasta"), total=args.total):
    mapper.add_genome(record.id, record.seq.encode())
    sizes[record.id] = len(record.seq)

print("Indexing...")
mapper.index()

ids = sorted(sizes)
ids_index = { name:i for i,name in enumerate(ids) }
ds = disjoint_set.DisjointSet({i:i for i in range(len(sizes)) })




print("Mapping...")

def mapping(record, pbar):
    hits = mapper.query_genome(record.seq.encode())
    pbar.update(1)
    return record.id, hits

with multiprocessing.pool.ThreadPool() as pool:
    with tqdm.tqdm(total=len(sizes)) as pbar:
        records = Bio.SeqIO.parse(args.input, "fasta")
        all_hits = pool.map(partial(mapping, pbar=pbar), records)
    for record_id, hits in all_hits:
        i = ids_index[record_id]
        for hit in filter(lambda hit: hit.identity >= args.cutoff, hits):
            j = ids_index[hit.name]
            ds.union(i, j)

contig_sets = list(ds.itersets())
print("Found", len(contig_sets), "independent sets")

print("Recovering largest sequence of each set")
representatives = set()
for contig_set in contig_sets:
    representative_i = max( contig_set, key=lambda i: sizes[ids[i]] )
    representatives.add(ids[representative_i])

assert len(representatives) == len(contig_sets)

# 6 genomes -> one file per genome
print("Sketching the 6 genomes")
mapper = pyfastani.Mapper()
for gbk in glob.glob(os.path.join(args.six, "*.gbk")):
    record = Bio.SeqIO.read(gbk, "genbank")
    mapper.add_genome(record.id, record.seq.encode())

# 9 genomes -> some fragment, manual handling needed
print("Sketching the 9 genomes")
for gbk in glob.glob(os.path.join(args.nine, "*.gbk")):
    base = os.path.basename(gbk)
    if not base.startswith(("Streptomyces_ghanaensis", "Streptomyces_sp_AA4", "Streptomyces_sp_C")):
        record = Bio.SeqIO.read(gbk, "genbank")
        mapper.add_genome(record.id, record.seq.encode())
# two contigs for S. ghanaensis
g1 = Bio.SeqIO.read(os.path.join(args.nine, "Streptomyces_ghanaensis_ATCC_14672_DS999641.1.gbk"), "genbank")
g2 = Bio.SeqIO.read(os.path.join(args.nine, "Streptomyces_ghanaensis_ATCC_14672_DS999642.1.gbk"), "genbank")
mapper.add_draft(g1.id, [g1.seq.encode(), g2.seq.encode()])
# two contigs for S. sp. AA4
a1 = Bio.SeqIO.read(os.path.join(args.nine, "Streptomyces_sp_AA4_GG657746.1.gbk"), "genbank")
a2 = Bio.SeqIO.read(os.path.join(args.nine, "Streptomyces_sp_AA4_GG657747.1.gbk"), "genbank")
mapper.add_draft(a1.id, [a1.seq.encode(), a2.seq.encode()])
# three contigs for S. sp. C
c1 = Bio.SeqIO.read(os.path.join(args.nine, "Streptomyces_sp_C_GG657750.1.gbk"), "genbank")
c2 = Bio.SeqIO.read(os.path.join(args.nine, "Streptomyces_sp_C_GG657751.1.gbk"), "genbank")
c3 = Bio.SeqIO.read(os.path.join(args.nine, "Streptomyces_sp_C_GG657752.1.gbk"), "genbank")
mapper.add_draft(c1.id, [c1.seq.encode(), c2.seq.encode(), c3.seq.encode()])

print("Indexing")
mapper.index()

print("Mapping...")
def mapping(record, pbar):
    hits = mapper.query_genome(record.seq.encode())
    pbar.update(1)
    return record.id, hits

with multiprocessing.pool.ThreadPool() as pool:
    with tqdm.tqdm(total=len(representatives)) as pbar:
        records = Bio.SeqIO.parse(args.input, "fasta")
        all_hits = pool.map(partial(mapping, pbar=pbar), filter(lambda r: r.id in representatives, records))

print("Blacklisting")
blacklist = {
    representative_id
    for representative_id, hits in all_hits
    if any(hit.identity > args.cutoff for hit in hits)
}
print("Blacklisted", len(blacklist), "contigs")

whitelist = representatives - blacklist
print("Writing list of", len(whitelist), "remaining contigs")
with open(args.list, "w") as f:
    f.writelines(f"{r}\n" for r in sorted(whitelist))

if args.output is not None:
    print("Writing output")
    with open(args.output, "w") as f:
        for i, record in enumerate(tqdm.tqdm(Bio.SeqIO.parse(args.input, "fasta"), total=len(sizes))):
            if record.id in whitelist:
                Bio.SeqIO.write(record, f, "fasta")
