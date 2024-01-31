#!/usr/bin/env python3
"""Count overlapping of VIP regions vs decoy/real selection regions."""
import sys
import numpy as np

in_region_filename = sys.argv[1]
in_decoy_filename = sys.argv[2]

# in_region_filename = 'data/cov_vip.regions'
# in_decoy_filename = 'data/ckb_real_decoy_sr_regions'


###############################################################################
# functions for overlapping count


def find_intersection_length(region1, region2):
    """Return the overlap length between two intervals."""
    r1l, r1r = region1[0], region1[1]
    r2l, r2r = region2[0], region2[1]
    if r1r < r2l:
        return 0
    if r2r < r1l:
        return 0
    ovlp_l, ovlp_r = max(r1l, r2l), min(r1r, r2r)
    return ovlp_r - ovlp_l + 1


def find_intersection(region1, region2):
    """Return the intersection left and right, if two intervals overlap, \
otherwise none."""
    r1l, r1r = region1[0], region1[1]
    r2l, r2r = region2[0], region2[1]
    if r1r < r2l:
        return None
    if r2r < r1l:
        return None
    return [max(r1l, r2l), min(r1r, r2r)]


def merge_wins(chr_wins):
    """Merge loci of a chromosome, given a list of [bp_start, bp_end]."""
    if not chr_wins:
        return []
    starts_and_ends = []
    for win in chr_wins:
        starts_and_ends.append((win[0], 'start'))
        starts_and_ends.append((win[1], 'end'))
    starts_and_ends = sorted(starts_and_ends, key=lambda x: x[0])
    merged_chr_wins = []
    stack = [starts_and_ends[0]]
    for i in range(1, len(starts_and_ends)):
        brk = starts_and_ends[i]
        if brk[1] == 'start':
            stack.append(brk)
        elif brk[1] == 'end':
            if stack:
                start = stack.pop()[0]
            if len(stack) == 0:
                end = brk[0]
                merged_chr_wins.append((start, end))
    return merged_chr_wins


###############################################################################
# read gene regions
print('read gene regions', file=sys.stderr)

inlines = open(in_region_filename).read().splitlines()
gene_regions = [[line.split()[0], int(line.split()[1]), int(line.split()[2])]
                for line in inlines]

# adding the + -10K
gene_regions = [[rg[0], max(1, rg[1]-10000), rg[2]+10000]
                for rg in gene_regions]

# they should be in the CKB genotyped region
inlines = open('data/ckb_regions').read().splitlines()
ckb_regions = [[line.split()[0], int(line.split()[1]), int(line.split()[2])]
               for line in inlines]

sum_gene_in_ckb = 0
sum_gene_bp = 0
gene_chr_regions = {}
for chr_int in range(1, 23):
    chr = f'{chr_int:d}'
    cg_regions = [[rg[1], rg[2]] for rg in gene_regions if rg[0] == chr]
    cc_regions = [[rg[1], rg[2]] for rg in ckb_regions if rg[0] == chr]
    all_chr_rgs = []
    for rg1 in cg_regions:
        is_gene_in_ckb = False
        for rg2 in cc_regions:
            ovlp_len = find_intersection_length(rg1, rg2)
            if not ovlp_len:
                continue
            is_gene_in_ckb = True
            sum_gene_bp += ovlp_len
            all_chr_rgs.append(find_intersection(rg1, rg2))
        if is_gene_in_ckb:
            sum_gene_in_ckb += 1

    all_chr_merged_rgs = merge_wins(all_chr_rgs)
    gene_chr_regions[chr] = []
    for rg in all_chr_merged_rgs:
        gene_chr_regions[chr].append([rg[0], rg[1]])


###############################################################################
# read the SR lines and count the overlaping

print('read SR regions and count overlap', file=sys.stderr)

sr_lines = open(in_decoy_filename).read().splitlines()

decoy_results = []  # ovlp_count, over_half_ovlp_count, ovlp_bp_count
for sr_line in sr_lines:
    sr_regions = sr_line.split(',')
    sr_regions = [[sr.split()[0], int(sr.split()[1]), int(sr.split()[2])]
                  for sr in sr_regions]
    sum_gene_in_sr = 0
    sum_gene_over_half_in_sr = 0
    sum_gene_in_sr_bp = 0
    for chr_int in range(1, 23):
        chr = f'{chr_int:d}'
        chr_sr_regions = [[rg[1], rg[2]] for rg in sr_regions if rg[0] == chr]
        for rg1 in gene_chr_regions[chr]:
            is_gene_in_sr = False
            is_gene_over_half_in_sr = False
            for rg2 in chr_sr_regions:
                ovlp_len = find_intersection_length(rg1, rg2)
                if ovlp_len == 0:
                    continue
                is_gene_in_sr = True
                sum_gene_in_sr_bp += ovlp_len
                rg1_len = rg1[1] - rg1[0] + 1
                rg2_len = rg2[1] - rg2[0] + 1
                if ovlp_len >= rg1_len*0.5 or ovlp_len >= rg2_len*0.5:
                    is_gene_over_half_in_sr = True
            if is_gene_in_sr:
                sum_gene_in_sr += 1
            if is_gene_over_half_in_sr:
                sum_gene_over_half_in_sr += 1
    decoy_results.append([sum_gene_in_sr, sum_gene_over_half_in_sr,
                          sum_gene_in_sr_bp])

###############################################################################
# sort, get the p values and print

real_ovlp_counts = decoy_results[0]
# gene_ovlp, geno_gt_half_ovlp, BP_ovlp

# the ovlp_count P value
ovlp_p = [1, 1, 1]
tile_0025 = [None, None, None]
tile_0500 = [None, None, None]
tile_0975 = [None, None, None]
ovlp_mean = [None, None, None]
ovlp_std = [None, None, None]

for i in range(3):
    ovlp_list = sorted([result[i] for result in decoy_results])[::-1]
    real_ovlp_count = real_ovlp_counts[i]
    p1 = ovlp_list.index(real_ovlp_count) + 1
    p1 += (ovlp_list.count(real_ovlp_count)-1) / 2
    p1 /= 10000.0
    ovlp_p[i] = p1
    tile_0025[i] = ovlp_list[250]
    tile_0500[i] = ovlp_list[5000]
    tile_0975[i] = ovlp_list[9750]
    ovlp_mean[i] = np.mean(ovlp_list)
    ovlp_std[i] = np.std(ovlp_list)

in_name = in_region_filename.replace('regions', '').replace('data/', '').\
    replace('_regions', '').replace('.', '').upper()
sr_name = in_decoy_filename.replace('data/', '').replace('real_decoy_', '').\
    replace('_regions', '').upper()

# output
out_str = f'{in_name} {sr_name} {sum_gene_in_ckb} {sum_gene_bp}'
for i in range(3):
    out_str += f' {ovlp_mean[i]:.3f} {ovlp_std[i]:.3f} \
{tile_0975[i]} {tile_0500[i]} {tile_0025[i]} \
{real_ovlp_counts[i]} {ovlp_p[i]}'

print(out_str.replace(' ', ','))

###############################################################################
