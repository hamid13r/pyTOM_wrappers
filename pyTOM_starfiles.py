#!/home/janlevinson/micromamba/envs/pytom_env/bin/python
# Author: Hamidreza Rahmani
# Grotjahn  Lab
# A quick script that works as a wrapper for pyTOM
#before running  this scrips first run: conda activate pytom_env
#This is part 1 of the wrapper that extract a high number of candidates with a low threshold
#After running this part the threshold needs to be determined with the help of 
#plotGaussianFit.py -f candidates/[filename].xml  -n 300 -p -1
#and then the second part of the wrapper to make the startfiles.
#Starfiles should be checked (preferably using cube) and change the threshold if necessary before extraction in Warp




import numpy as np
import tqdm
import click
import glob
import os
from pathlib import Path
import xml.etree.ElementTree as ET
import sys
import pandas as pd  

@click.command()
@click.option('--candidates_folder', help='Path to the folder for candidate xml', required=False, default="candidates/",show_default=True)
@click.option('--threshold_folder', help='Path to the folder for candidate xml', required=False, default="candidates_threshold/",show_default=True)
@click.option('--star_folder', help='Path to the folder for candidate xml', required=False, default="starfiles/",show_default=True)
@click.option('--plots_folder', help='Path to the folder for histograms plots', required=False, default="plots/",show_default=True)
@click.option('--ref', help='Reference Name', required=True,show_default=True)
@click.option('--threshold', help='Minimum score for candidates to be extracted', required=False, type=float, default=0.05,show_default=True)
@click.option('--auto_threshold', help='Choose if you want threshold to be calculated by pyTOM.', required=False, is_flag=True, default=False,show_default=True)
@click.option('--apix', help='Pixel size', required=True)
@click.option('--z_list', help='A file with the list of min/max Z coordinates for each tomogram, tab-serparated, tomogram zmin zmax', default="zlist.txt", required=False,show_default=True)



def main(candidates_folder, threshold_folder, plots_folder, ref, auto_threshold, star_folder, threshold, apix, z_list):

	#os.system("conda activate pytom_env")
	if not os.path.isdir(candidates_folder):
		print("Candidates directory not found, STOPPING")
		sys.exit(0)
	if not os.path.isdir(star_folder):
		print("Candidates direcory not found, making directory %s" %(star_folder))
		os.mkdir(star_folder)
	if not os.path.isdir(threshold_folder):
		print("Candidates direcory not found, making directory %s" %(threshold_folder))
		os.mkdir(threshold_folder)
	if not os.path.isdir(plots_folder):
		print("Candidates direcory not found, making directory %s" %(plots_folder))
		os.mkdir(plots_folder)
	if os.path.isfile(z_list):
		zdf = pd.read_csv(z_list, delimiter="\t", header=None, names=['tomo','z_min','z_max'])

	refname = Path(ref).stem
	for xmlfile in glob.glob(f"{candidates_folder}/*{refname}.xml"):
		xmlname = Path(xmlfile).stem

		mydoc = ET.parse(xmlfile)
		root = mydoc.getroot()
		
		if auto_threshold:
			threshold = os.popen(f"plotGaussianFit.py -f {xmlfile} -n 200 -p -1 -o {plots_folder}/{xmlname}.png | grep \"optimal correlation coefficient threshold\" ").read()
			threshold = threshold.split(" ")[-1]
			threshold = float(threshold[:-2])
		
		tomoname = xmlname.replace("_"+refname,'') + ".mrc"
		if os.path.isfile(z_list):
			z_min = zdf.loc[zdf['tomo'] == tomoname]['z_min'].values[0]
			z_max = zdf.loc[zdf['tomo'] == tomoname]['z_max'].values[0]
		else:
			print("z_list not used.")
			z_min = -999
			z_max = 9999


		for child in root.findall('Particle'):
			if float(child[5].attrib['Value']) < threshold:
				root.remove(child)	
			elif int(child[2].attrib['Z']) < z_min:
				root.remove(child)
			elif int(child[2].attrib['Z']) > z_max:
				root.remove(child)
			
		

		out_tree = ET.ElementTree(root)
		with open("%s/%s.xml" %(threshold_folder,xmlname),'wb') as f:
			out_tree.write(f)

		
		

		convert_command=f"convert.py -f {threshold_folder}/{xmlname}.xml -t {star_folder} --outname {xmlname}.star -o star --pixelSize {apix} --binPyTom 1 --binWarpM 1"
		os.popen(convert_command)
		particle_numder = len(root)
		print(f"{threshold_folder}/{xmlname}.xml converted to star file at threshold of {threshold}, {particle_numder} particles located.")
		
			
if __name__ == '__main__':
    main()