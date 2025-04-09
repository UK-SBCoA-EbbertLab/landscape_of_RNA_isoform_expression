#!/usr/bin/env python

print("inside python script", flush=True)
import sys
import biolib
biolib.login()

print("imports done and logged in", flush=True)

biolib.utils.STREAM_STDOUT = True
print("loading tool", flush=True)
deeptmhmm = biolib.load('DTU/DeepTMHMM:1.0.24')
print("about to run", flush=True)
deeptmhmm_job = deeptmhmm.cli(args=f"--fasta {sys.argv[1]}")
print("saving files", flush=True)
deeptmhmm_job.save_files('result')
print("done in python", flush=True)
