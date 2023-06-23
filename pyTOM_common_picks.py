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
import pandas as pd
import glob
import starfile
import click
from scipy.spatial import distance

@click.command()
@click.option('--input_dir', prompt='Input directory', help='The directory of input files.')
@click.option('--ref1', prompt='Reference A name', help='The name of reference A.')
@click.option('--ref2', prompt='Reference B name', help='The name of reference B.')
@click.option('--threshold', prompt='Threshold', help='The threshold for Euclidean distance.', type=float, default=5)
def process_files(input_dir, ref1, ref2, threshold):
    """
    This function processes .star files in a given directory, comparing 3D coordinates between two references (ref1 and ref2) in each file.
    It calculates the Euclidean distance between the 3D coordinates and identifies coordinates that are closer than a given threshold.
    The coordinates that meet this condition are saved to new .star files.

    Parameters:
    input_dir (str): The directory of input files.
    ref1 (str): The name of reference A.
    ref2 (str): The name of reference B.
    threshold (float): The threshold for Euclidean distance.

    Example:
    To run this function, use the command:
    python /home/janlevinson/bin/pyTOM_common_picks.py --input_dir /path/to/input/files --ref1 referenceA --ref2 referenceB --threshold 0.5
    """
    ref1 = ref1.replace(".mrc","")
    ref2 = ref2.replace(".mrc","")
    # Get list of all .star files
    star_filesA = glob.glob(f'{input_dir}/*_{ref1}.star')
    #star_filesB = glob.glob(f'{input_dir}/*_{ref2}.star')

    # Read each file into a pandas dataframe
    for sfile in star_filesA:
        dfA = starfile.read(sfile)
        
        dfB = starfile.read(sfile.replace(f'_{ref1}.star', f'_{ref2}.star'))

        # Compare 3D coordinates in each dataframe
        # Calculate Euclidean distance between 3D coordinates
        dist = distance.cdist(dfA[['rlnCoordinateX', 'rlnCoordinateY', 'rlnCoordinateZ']], dfB[['rlnCoordinateX', 'rlnCoordinateY', 'rlnCoordinateZ']], 'euclidean')
        # Get indices of coordinates closer than threshold in dfA
        indicesA = np.where(dist < threshold)
        # Get indices of coordinates closer than threshold in dfB
        indicesB = np.where(dist.T < threshold)
        # Get coordinates in dfA that are found in dfB
        dfA_inB = dfA.iloc[indicesA[0]]
        # Get coordinates in dfB that are found in dfA
        dfB_inA = dfB.iloc[indicesB[0]]
        # Save coordinates to new .star files
        starfile.write(dfA_inB, sfile.replace(f'_{ref1}.star', f'_{ref1}_common.star'), overwrite=True)
        starfile.write(dfB_inA, sfile.replace(f'_{ref1}.star', f'_{ref2}_common.star'), overwrite=True)
if __name__ == '__main__':
    process_files()
