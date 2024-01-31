#!/usr/bin/env python3
#$ -cwd -V
#$ -N loading_parse -j y
#$ -q short.qc
#$ -l h_vmem=13G

import sys
import numpy as np
import scipy
import scipy.stats
import subprocess

################################################################################
# HMM parsing of SNP loadings                                                  #
# during the iterative high-LD region removal                                  #
################################################################################

################################################################################
# parameters

WING_SIZE = 500000
# adding this much to to the starts and ends of chromosomes
# for the calculation of cumulative recombination rates.

CRR_SCALE = 1E+7
# transition probability between states
# t_p = crr/CRR_SCALE
# transp = np.array([[1-t_p, t_p], [t_p, 1-t_p]])
# 7 wins in PC1-8 with 1E+8 (iter 02)
# 7 wins in PC1-8 with 1E+7 (iter 02) as well

# emission probability calculation:

LOADING_SCALE = 5
# e1 = scipy.stats.chi2.cdf(loading*loading,df=1,scale=LOADING_SCALE)

################################################################################
# arguments

iter_int = int(sys.argv[1])
# iter_int = 4
iter = f'{iter_int:02d}'

# parse how many PCs
# MAX_PC = int(sys.argv[2])
MAX_PC = 20

################################################################################
# get the cumulative recombination rates

rcomb_lines = open('data/rcomb_sum.txt').read().splitlines()

inlines=open('snp_loading_in/iter_'+iter+'.bim').read().splitlines()
snps     = [line.split()[1] for line in inlines ]
snp2indx = {snps[i]:i for i in range(len(snps))}
snp2chr  = {line.split()[1]:line.split()[0] for line in inlines}
snp2bp   = {line.split()[1]:int(line.split()[3]) for line in inlines}

chr_trans ={}

for chr_int in range(1,23):
    chr = f'{chr_int:d}'
    chr_snps = [snp for snp in snps if snp2chr[snp]==chr ]
    chr_bps  = [snp2bp[snp] for snp in chr_snps ]
    chr_bp2snp = {chr_bps[i]:chr_snps[i] for i in range(len(chr_bps))}
    #
    start_bp = chr_bps[0] - WING_SIZE
    if start_bp < 1 : start_bp = 1
    end_bp = chr_bps[-1] + WING_SIZE
    #
    chr_rcomb_lines = [line for line in rcomb_lines  if line.split()[0] == chr]
    chr_bp2srb ={int(line.split()[2]) : float(line.split()[1]) \
        for line in rcomb_lines if line.split()[0] == chr}
    chr_bp2srb = {bp: chr_bp2srb[bp] for bp in chr_bp2srb.keys() \
        if bp > start_bp and bp < end_bp }
    #
    for bp in chr_bps:
        if bp not in chr_bp2srb.keys():
            chr_bp2srb[bp] = 0.0
    #
    chr_bp_set = set(chr_bps)
    bp_index={}
    cur_index=0
    for bp in chr_bp2srb.keys():
        bp_index[bp]=cur_index
        if bp in chr_bp_set:
            cur_index += 1
    ####################################
    #cumulative recombination rates
    chr_bps.append(end_bp)
    crrs=[0.0]*len(chr_bps)
    for bp in chr_bp2srb.keys():
        index = bp_index[bp]
        crrs[index] += chr_bp2srb[bp] #
        if bp == chr_bps[index] and bp < chr_bps[-2]:
            crrs[index]   -= chr_bp2srb[bp] * 0.5
            crrs[index+1] += chr_bp2srb[bp] * 0.5
    # crr is between 0 and 135552, mean 183, sdv 450 (iter 2)
    chr_trans[chr] = []
    for crr in crrs:
        t_p = crr/CRR_SCALE # transition probability between stats
        transp = np.array([[1-t_p,t_p],[t_p,1-t_p]])
        chr_trans[chr].append(transp)

# 11 seconds so far
################################################################################
# HMM parsing, calculate marginal likelihood

inlines=open("snp_loading_in/iter_"+iter+".eigenvec.var").read().splitlines()

# default is the plink gcta output
is_flashpca_output = False
if inlines[0].startswith('SNP\tRefAllele\tV1'):
    is_flashpca_output=True
    inlines.pop(0)

head_col = ['SNP','CHR','BP']
for pc in range(MAX_PC):
    head_col.append('PC'+str(pc+1)+'_loading')
    head_col.append('PC'+str(pc+1)+'_high_ld_p')

out_loadings=[]
out_high_ld_ps=[]
for snp in snps:
    out_loadings.append([])
    out_high_ld_ps.append([])

for pc in range(MAX_PC):
    pc_loadings = []
    if is_flashpca_output:
        pc_loadings = [float(line.split()[2+pc]) for line in inlines]
    else:
        pc_loadings = [float(line.split()[4+pc]) for line in inlines]
    ####################################
    # transfer loading into standard normal
    pc_loadings = scipy.stats.zscore(pc_loadings)
    #
    for chr_int in range(1,23):
        chr = f'{chr_int:d}'
        chr_snps     = [ snp for snp in snps if snp2chr[snp]== chr ]
        chr_loadings = [ pc_loadings[snp2indx[snp]] for snp in chr_snps ]
        chr_loadings = np.array(chr_loadings)
        ################################
        # emission probabilities
        e1 = scipy.stats.chi2.cdf(chr_loadings*chr_loadings, \
                                  df = 1, scale = LOADING_SCALE)
        e0 = 1 - e1
        transps = chr_trans[chr]
        ################################
        # forward
        fw_ps = [[0.99,0.01]]
        for i in range(len(e0)):
            fw_p0 =  fw_ps[i][0] * transps[i][0][0] \
                   + fw_ps[i][1] * transps[i][1][0]
            fw_p0 *= e0[i]
            fw_p1 =  fw_ps[i][0] * transps[i][0][1] \
                   + fw_ps[i][1] * transps[i][1][1]
            fw_p1 *= e1[i]
            fw_p0, fw_p1 = fw_p0/(fw_p0+fw_p1) , fw_p1/(fw_p0+fw_p1)
            fw_ps.append([fw_p0, fw_p1])
        ################################
        # backward
        bw_ps = [[0.99,0.01]]
        for i in range(len(e0)):
            min_i = len(e0) - i - 1
            bw_p0 =    bw_ps[i][0] * transps[min_i+1][0][0]  \
                     + bw_ps[i][1] * transps[min_i+1][1][0]
            bw_p0 *= e0[min_i]
            bw_p1 =    bw_ps[i][0] * transps[min_i+1][0][1]  \
                     + bw_ps[i][1] * transps[min_i+1][1][1]
            bw_p1 *= e1[min_i]
            bw_p0, bw_p1 = bw_p0/(bw_p0+bw_p1) , bw_p1/(bw_p0+bw_p1)
            bw_ps.append([bw_p0,bw_p1])
        ################################
        # marginal
        marg_p1=[]
        for i in range(len(e0)):
            m_p0, m_p1 = fw_ps[i+1][0],fw_ps[i+1][1]
            m_p0 *= bw_ps[len(e0)-i-1][0]
            m_p1 *= bw_ps[len(e0)-i-1][1]
            m_p0, m_p1 = np.sqrt(m_p0), np.sqrt(m_p1)
            m_p0, m_p1 = m_p0 / (m_p0 + m_p1) , m_p1 / (m_p0 + m_p1)
            marg_p1.append(m_p1)
        ################################
        # output
        for i in range(len(e0)):
            indx = snp2indx[chr_snps[i]]
            out_loadings  [indx].append(np.abs(chr_loadings[i]))
            out_high_ld_ps[indx].append(marg_p1[i])


ofile_name='snp_loading_out/iter_'+iter+'_loadings.txt'
ofile=open(ofile_name,'w')
print(' '.join(head_col),file=ofile)
for snp in snps:
    indx = snp2indx[snp]
    out_col = [snp,snp2chr[snp],str(snp2bp[snp])]
    for pc in range(MAX_PC):
        out_col.append('{:.5g}'.format(out_loadings  [indx][pc]))
        out_col.append('{:.5g}'.format(out_high_ld_ps[indx][pc]))
    print(' '.join(out_col),file=ofile)
ofile.close()

# 48 seconds so far
# exit(0)

################################################################################
# plot loading

plot_cmd = 'iter_loading_plot.R {}'.format(iter_int)
subprocess.run(plot_cmd,shell=True)

################################################################################
