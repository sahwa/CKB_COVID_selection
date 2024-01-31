#!/usr/bin/env python3
"""Use bedtools to generate randomly permuted SR regions."""

import subprocess
import sys

in_sr_region_filename = sys.argv[1]
chr = sys.argv[2]
out_decoy_file_root = sys.argv[3]
# make_decoy.py data/ukb_sr_regions 1 decoy/ukb

inlines = open('data/ckb_regions').read().splitlines()
ckb_regions = [[int(line.split()[1]), int(line.split()[2])]
               for line in inlines if line.split()[0] == chr]

inlines = open(in_sr_region_filename).read().splitlines()
sr_regions = [[int(line.split()[1]), int(line.split()[2]), line.split()[3]]
              for line in inlines if line.split()[0] == chr]

if not sr_regions:
    exit(0)

in_sr_bed_file = open(f'in_chr{chr}.bed', 'w')
for rg in sr_regions:
    print(f'chr{chr}\t{rg[0]-1}\t{rg[1]}\t{rg[2]}\t0\t+', file=in_sr_bed_file)
in_sr_bed_file.close()

in_genome_file = open(f'in_chr{chr}.genome', 'w')
print(f'chr{chr}\t{ckb_regions[-1][1]}', file=in_genome_file)
in_genome_file.close()

in_exclu_file = open(f'in_chr{chr}_exclude.bed', 'w')
if ckb_regions[0][0] != 1:
    print(f'chr{chr}\t0\t{ckb_regions[0][0]-1}', file=in_exclu_file)

for i, rg in enumerate(ckb_regions):
    if i == 0:
        continue
    print(f'chr{chr}\t{ckb_regions[i-1][1]}\t{rg[0]-1}', file=in_exclu_file)
in_exclu_file.close()

for ll in range(10000):
    ofile_name = f'{out_decoy_file_root}_{chr}_{ll+1:06d}'
    cmd = f'bedtools shuffle -i in_chr{chr}.bed \
        -excl in_chr{chr}_exclude.bed \
        -g in_chr{chr}.genome \
        -seed {ll+1} \
         -noOverlapping > {ofile_name}'
    print(cmd, file=sys.stderr)
    subprocess.call(cmd, shell=True)
