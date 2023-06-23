#!/home/janlevinson/micromamba/envs/pytom_env/bin/python
# Author: Hamidreza Rahmani
# Grotjahn  Lab
# A quick script that works as a wrapper for pyTOM
#before running  this scrips first run: conda activate pytom_env
#This is part 1 of the wrapper that extract a high number of candidates with a low threshold
#After running this part the threshold needs to be determined, we normally do it by changing the threshold in the wrapper and
#looking at the star file in Cube. 
#and then the second part of the wrapper to make the startfiles.
#Starfiles should be checked (preferably using cube) and change the threshold if necessary before extraction in Warp




import numpy as np
import tqdm
import click
import glob
import os
from pathlib import Path
import sys
import subprocess
import pandas as pd
from multiprocessing import Pool


     

@click.command()
@click.option('--tomo_folder', help='Path to the folder containing tomograms', required=True, default="tomo/",show_default=True)
@click.option('--tomomask_folder', help='Path to the folder containing tomograms masks', required=False, default=None,show_default=True)
@click.option('--cc_folder', help='Path to the folder containing cross-corelation maps', required=True, default="cc/",show_default=True)
@click.option('--job_folder', help='Path to the folder containing jobfiles', required=True, default="jobfiles/",show_default=True)
@click.option('--candidates_folder', help='Path to the folder for candidate xml', required=True, default="candidates/",show_default=True)
@click.option('--candidates_number', help='Number of candidates to be extracted', required=True, type=int, default=5000,show_default=True)
@click.option('--command_file_prefix', help='File prefix to write commands in.', required=False, default="commands",show_default=True)
@click.option('--extraction_threads', help='Number of parallel extractions to run', required=True, type=int, default=48,show_default=True)
@click.option('--min_score', help='Minimum score for candidates to be extracted', required=True, type=float, default=0.05,show_default=True)
@click.option('--ref', help='Path to reference file (BonW, same apix)', required=True,show_default=True)
@click.option('--mask', help='Path to the mask applied to reference', required=True,show_default=True)
@click.option('--apix', help='Pixel size', required=True,show_default=True)
@click.option('--z_min', help='Minimum Z to find particles in tomograms', type=int, required=True,show_default=True)
@click.option('--z_max', help='Maximum Z to find particles in tomograms', type=int, required=True,show_default=True)
@click.option('--wedge1', help='Wedge angle (90-maxtilt)', required=True, type=int, default=30,show_default=True)
@click.option('--wedge2', help='Wedge angle (90-maxtilt)', required=True, type=int, default=30,show_default=True)
@click.option('--min_distance', help='Minimum distance between candidates', required=True, type=int, default=10,show_default=True)
@click.option('--margin', help='Margins from edges of tomogram for candidates', required=True, type=int, default=10,show_default=True)
@click.option('--anglist', help='Angle list file', required=True, default="angles_19.95_1944.em",show_default=True)
@click.option('--skip_cc', help='Skip cross-correlation calculation and go to extracting candidates.', required=False, is_flag=True, default=False,show_default=True)
@click.option('--dryrun', help='Just make command files and do not run anything...', required=False, is_flag=True, default=False,show_default=True)
def main(tomo_folder, tomomask_folder, cc_folder, job_folder, candidates_folder, candidates_number, command_file_prefix, extraction_threads, min_score, ref, mask, apix, wedge1, wedge2, z_min, z_max, anglist, min_distance, margin, skip_cc, dryrun):
	"""Example: pyTOM_candidates.py --tomo_folder tomo/ --tomomask_folder tomomask/ --cc_folder cc/ --job_folder jobfiles/ --candidates_folder candidates/ --candidates_number 5000 --command_file_prefix extraction_commands.txt
     --extraction_threads 48 --min_score 0.05 --ref 50Sref10.mrc --mask 50Smask10.mrc --apix 9.98 --z_min 0 --z_max 266 --wedge1 30 --wedge2 30 --min_distance 10 --margin 10 --anglist angles_19.95_1944.em --skip_cc"""
	if not os.path.isdir(tomo_folder):
		print("Tomogram directory not found, STOPPING")
		sys.exit(0)
	if not (tomomask_folder is None):
		if not os.path.isdir(tomomask_folder):
			print("Tomogram MASKS directory not found, STOPPING")
			sys.exit(0)
	if not os.path.isdir(cc_folder):
		print("CC directory not found, making directory %s" %(cc_folder))
		os.mkdir(cc_folder)
	if not os.path.isdir(job_folder):
		print("Jobs directory not found, making directory %s" %(job_folder))
		os.mkdir(job_folder)
	if not os.path.isdir(candidates_folder):
		print("Candidates directory not found, making directory %s" %(candidates_folder))
		os.mkdir(candidates_folder)
			
	counter = 0
	tomolist = glob.glob("%s/*.mrc" %(tomo_folder))
	tomonum = len(tomolist)
	#if os.path.isfile(z_list):
	#	for tomofile in tomolist:
	#		if not (zdf['tomo'].eq(os.path.basename(tomofile))).any():
	#			print(f"Z values for {os.path.basename(tomofile)} not found, make sure you have the file formatted correctly... exisitng!")
	#			sys.exit(0)
	job_commands=[]
	for tomofile in tomolist:
		tomoname = Path(tomofile).stem
		refname = Path(ref).stem
		tomo_cc_folder = cc_folder + tomoname
		if not os.path.isdir(tomo_cc_folder):
			print("Making %s to save score and angle maps..." %(tomo_cc_folder))
			os.mkdir(tomo_cc_folder)

		#if os.path.isfile(z_list):
		#	z_min = zdf.loc[zdf['tomo'] == os.path.basename(tomofile)]['z_min'].values[0]
		#	z_max = zdf.loc[zdf['tomo'] == os.path.basename(tomofile)]['z_max'].values[0]
			#print("Z range read from file to be {z_min}-{z_max}, replacing the default values!")
		job_command=f"localizationJob.py -v {tomofile} -r {ref} -m {mask} --wedge1 {wedge1} --wedge2 {wedge2} -a {anglist}  -d {cc_folder}/{tomoname} --splitX 8 --splitY 12  --splitZ 4 -j {job_folder}/{tomoname}_{refname}.xml --zstart {z_min} --zend {z_min+z_max}"
		print(f"Making Job File: {job_folder}/{tomoname}_{refname}.xml")
		job_commands.append(job_command)
		if not dryrun:
			subprocess.getoutput(job_command)
	print("Job files are created.")
	
	with open(f"{command_file_prefix}_{refname}_jobs.txt", 'w') as f:
		for command in job_commands:
			f.write("%s\n" % command)       
	print(f"\n\nJobfile creating commands are written to {command_file_prefix}_{refname}_jobs.txt")
	
	cc_commands=[]
	if not skip_cc:		
		for tomofile in tomolist:
			tomoname = Path(tomofile).stem
			refname = Path(ref).stem
			tomo_cc_folder = cc_folder + tomoname
			counter = counter + 1
			print(f"\n\nCross-Correlating TOMOGRAM {counter} OF {tomonum} USING GPU")
			cc_command=f"localization.py -j {job_folder}/{tomoname}_{refname}.xml -g 0"
			print(f"Calculating Cross-Correlations for {job_folder}/{tomoname}_{refname}.xml")
			if not dryrun:
				os.system(cc_command)
			cc_commands.append(cc_command)
		print(f"\n\nCross-Correlating is DONE, GPUs wont be used anymore!")
		with open(f"{command_file_prefix}_{refname}_cc.txt", 'w') as f:
			for command in cc_commands:
				f.write("%s\n" % command)
		print(f"\n\nCross-Correlation commands are written to {command_file_prefix}_{refname}_cc.txt")
	else:
		print(f"\n\nCross-Correlating was Skipped, no GPU usage...")	
		
	
	commands=[]
	counter = 0
	for tomofile in tomolist:
		tomoname = Path(tomofile).stem
		refname = Path(ref).stem
		tomo_cc_folder = cc_folder + tomoname

		counter = counter + 1
		if (tomomask_folder is None):
			extract_command=f"extractCandidates.py -j {job_folder}/{tomoname}_{refname}.xml -r {cc_folder}/{tomoname}/scores_{refname}.em -o {cc_folder}/{tomoname}/angles_{refname}.em -n {candidates_number} -s {min_distance} -p {candidates_folder}/{tomoname}_{refname}.xml -v {min_score} --margin {margin}" 
			commands.append(extract_command)
			#print(extract_command)
		else:
			extract_command=f"extractCandidates.py -j {job_folder}/{tomoname}_{refname}.xml --tomogram-mask {tomomask_folder}/{tomoname}.mrc -r {cc_folder}/{tomoname}/scores_{refname}.em -o {cc_folder}/{tomoname}/angles_{refname}.em -n {candidates_number} -s {min_distance} -p {candidates_folder}/{tomoname}_{refname}.xml -v {min_score} --margin {margin}" 
			commands.append(extract_command)
			#print(extract_command)
	with open(f"{command_file_prefix}_{refname}_extraction.txt", 'w') as f:
		for command in commands:
			f.write("%s\n" % command)       
	print(f"\n\nExtraction commands are written to {command_file_prefix}_extraction.txt")
	print(f"\n\nRunning extraction commands in parallel with {extraction_threads} threads, this will take a couple of hours for every batch of {extraction_threads} tomograms.")
	if not dryrun:
		os.system(f"cat {command_file_prefix}_{refname}_extraction.txt | parallel -j {extraction_threads} --jobs {extraction_threads}")
	print(f"\n\nExtraction is DONE, check {candidates_folder} for candidates.")
	print(f"\n\nYou can run the extraction commands in {command_file_prefix}_{refname} by running the following command:\n\nparallel -j {extraction_threads} < {command_file_prefix}_{refname}\n\n")
	

if __name__ == '__main__':
    main()
