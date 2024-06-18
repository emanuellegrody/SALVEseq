'''
Note about the script:
This is a python script to submit a sample to starcode to match for optimal Levenstein distance and merge similar sequences to ID all viral sequences
Output files:
	1. 3 files are created per sample where merged targets of varying lengths (18, 24, 30) are written

command to run this script: python3 starcodeRun.py
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''
from argparse import ArgumentParser
import os
import subprocess


starcodeLocation = os.path.join("/projects/b1042/GoyalLab/egrody/packages/starcode","starcode")

os.chdir("..")

lengths = [18, 24, 30]
samples = ["invitro"]
threads = "4"
dist = "8"

for i in lengths:
	for sample in samples:
		samplePath = os.path.join('/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/analysis/', sample, '{}_shavedReadsList_{}.txt'.format(sample, i))
		outfilePath = os.path.join('/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/analysis/starcode', sample, '{}_Target{}_d8.txt'.format(sample, i))
		starcodeCommand = [starcodeLocation, '-d', dist, '-t', threads, '-i', samplePath, "-s", '-o', outfilePath]
		subprocess.run(starcodeCommand)

