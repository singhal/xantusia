import re
import argparse
import os
import shutil
import pandas as pd
import numpy as np
import random

parser = argparse.ArgumentParser(description = "make SNP input file for moments",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# then get group 1 & 2
parser.add_argument("--sp1", type = str, default = None, help = "species 1 to compare")
parser.add_argument("--sp2", type = str, default = None, help = "species 2 to compare")
args = parser.parse_args()

vcffile = '/media/babs/brains3/Xantusia/data/no_out.vcf'
indfile = '/media/babs/brains3/Xantusia/data/xantusia_samples_v7.csv'
d = pd.read_csv(indfile)

sp1 = args.sp1
sp2 = args.sp2

indsp = {sp1: d.loc[d.OTU == sp1, "sample_fix"].tolist(),
         sp2: d.loc[d.OTU == sp2, "sample_fix"].tolist()}

snps = {}
f = open(vcffile, 'r')
for l in f:
    if re.search('#CHROM', l):
        inds = re.split("\t", l.rstrip())[9:]
    else:
        if not re.search("^#", l):
            d = re.split('\t', l.strip())
            alleles = [d[3]] + re.split(',', d[4])

            genos = [re.search('^(\S\S\S)', x).group(1) for x in d[9:]]
            genos = [re.split('/', x) for x in genos]
            
            genos = dict(zip(inds, genos))

            snp = {}
            for sp in indsp:
                snp[sp] = []
                for ind in indsp[sp]:
                    gen = genos[ind]
                    if gen[0] != '.':
                        gen = [alleles[int(x)] for x in gen]
                        snp[sp] += gen

            complete1 = len(snp[sp1]) / (2 * len(indsp[sp1]))
            complete2 = len(snp[sp2]) / (2 * len(indsp[sp2]))
            # check not too missing
            if complete1 > 0.5 and complete2 > 0.5:
                a = list(set(snp[sp1] + snp[sp2]))
                # check only two alleles
                # check actually variable
                if len(a) == 2:
                    # keep it!
                    if d[0] not in snps:
                        snps[d[0]] = {}
                    snps[d[0]][d[1]] = snp
f.close()

            
outf = 'data/moments/%s.%s.snps.txt' % (sp1, sp2)
of = open(outf, 'w')
of.write('Ingroup\tOutgroup\tAllele1\t%s\t%s\tAllele2\t%s\t%s\tGene\tPosition\n' % (sp1, sp2, sp1, sp2))
# we can't polarize
for c in snps:
    # pick a random snp to avoid LD
    pos = random.choice(list(snps[c].keys()))
    ingroup = snps[c][pos][sp1] + snps[c][pos][sp2]
    alleles = list(set(ingroup))

    inallele = alleles[0]
    outallele = alleles[1]
    outd = [None] * 10
    outd[0] = '-%s-' % inallele
    outd[1] = '-%s-' % outallele
    outd[2] = inallele
    outd[3] = str(snps[c][pos][sp1].count(inallele))
    outd[4] = str(snps[c][pos][sp2].count(inallele))
    outd[5] = outallele
    outd[6] = str(snps[c][pos][sp1].count(outallele))
    outd[7] = str(snps[c][pos][sp2].count(outallele))
    outd[8] = c
    outd[9] = str(pos)
    of.write('\t'.join(outd) + '\n')
of.close()

ct1 = len(indsp[sp1]) * 2
ct2 = len(indsp[sp2]) * 2

if ct1 < 3:
    ct1 = 4
if ct2 < 3:
    ct2 = 4

outsfs = '/media/babs/brains3/Xantusia/moments/input_files/%s.%s.sfs' % (sps['sp1'], sps['sp2'])
basesc = '/media/babs/brains3/Xantusia/scripts/run_moment_base.py'
with open(basesc, 'r') as file:
    filedata = file.read()
filedata = filedata.replace("snps = \"XXX\"", "snps = \"%s\"" % outf)
filedata = filedata.replace("fs.to_file(\"XXX\"", "fs.to_file(\"%s\"" % outsfs)
filedata = filedata.replace("pop_ids=[\"XXX\", \"XXX\"]", "pop_ids=[\"%s\", \"%s\"]" % (sps['sp1'], sps['sp2']))
filedata = filedata.replace("proj = [XXX, XXX]", "proj = [%s, %s]" % (ct1, ct2))
outdir = '/media/babs/brains3/Xantusia/moments/run_scripts/%s.%s' % (sps['sp1'], sps['sp2'])
if not os.path.isdir(outdir):
    os.mkdir(outdir)
newfile = '%s/run_moment.%s.%s.py' % (outdir, sps['sp1'], sps['sp2'])
shutil.copy('/media/babs/brains3/spheno/demography/moments/Models_2D.py', outdir)
shutil.copy('/media/babs/brains3/spheno/demography/moments/Optimize_Functions.py', outdir)
with open(newfile, 'w') as file:
  file.write(filedata)