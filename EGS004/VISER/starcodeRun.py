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


starcodeLocation = os.path.join("/home/egy2296/packages/starcode/","starcode")

os.chdir("..")

lengths = [18, 24, 30]
samples = ["VISER_D13", "VISER_D195", "VISER_D83"]
threads = "8"
dist = "8"

for i in lengths:
	for sample in samples:
		samplePath = os.path.join('/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS004/starcode/{}shavedReadsList_{}.txt'.format(sample, i))
		outfilePath = os.path.join('/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS004/starcode/{}_Target{}_d8.txt'.format(sample, i))
		starcodeCommand = [starcodeLocation, '-d', dist, '-t', threads, '-i', samplePath, "-s", '-o', outfilePath]
		subprocess.run(starcodeCommand)

