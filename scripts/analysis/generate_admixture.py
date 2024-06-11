import re
import os
import random
import subprocess
import glob
import argparse

parser = argparse.ArgumentParser(description="create pop gen runs")
parser.add_argument('--indir',
			default = None,
			help = "input dir")
parser.add_argument('--prefix',
			default = None,
			help = "prefix")
args = parser.parse_args()

# do admixture
files = glob.glob("%s/%s*ped" % (args.indir, args.prefix))
for file in files:
	file = re.sub('.ped', '', file)
	subprocess.call("~/bin/plink_1.9/plink --file %s --make-bed --noweb --out %s" % (file, file), shell = True)
	newfile = file + '.bed'
	for k in range(1, 13):
		base = re.sub('.ped', '', file)
		base = re.sub(".*/", "", base)
		print(base, k)
		subprocess.call("~/bin/admixture --cv %s %s > %s.cv.%s" % (newfile, k, base, k), shell = True)
		# os.rename("%s.%s.P" % (base, k), "%s.%s.%s.P" % (base, k, x))
		# os.rename("%s.%s.Q" % (base, k), "%s.%s.%s.Q" % (base, k, x))