#!/usr/bin/Python
# this script makes Dx-Relapse and Same-Phase folders and organizes files into respective folders.

import pandas as pd
import os
import subprocess
import glob 

files = glob.glob("*.txt")


dx_rel = []
same_phase = []

with open("all_models.txt") as f1:
	for line in f1:
		l = line.split("\t")
		for x in files:
			if(l[0] in x and l[3] in x):
				if(x in dx_rel):
					continue
				else:
					dx_rel.append(x)

for a in dx_rel:
	subprocess.call(["mv", a, "./Dx-Relapse"])
	


with open("same_phase_all_models.txt") as f2:
	for line in f2:
		l = line.split("\t")
		for x in files:
			if(l[0] in x and l[3] in x):
				if(x in same_phase):
					continue
				else:
					same_phase.append(x)


for b in same_phase:
	subprocess.call(["mv", b, "./Same-Phase"])




