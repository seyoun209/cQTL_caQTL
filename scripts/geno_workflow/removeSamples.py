#!/usr/bin/env python3

import sys
import pandas as pd
import difflib
import os

# Read in samples to remove
removeSamples = pd.read_csv(sys.argv[1], sep=" ")
removePLINK = []

for row in removeSamples.itertuples():
    # Grab fam file for batch
    fam_file = 'geno_output/binary/' + row.Batch + '.fam'
    
    # Check if file exists
    if not os.path.exists(fam_file):
        print(f"Warning: {fam_file} not found, skipping {row.Donor}")
        continue
        
    famFile = pd.read_csv(fam_file, sep=" ", header=None)
    
    # Match donor to 2nd column
    matches = difflib.get_close_matches(row.Donor, famFile.iloc[:,1], 1, cutoff=0.4)
    
    if matches:
        donorID = matches[0]
        # Select FID and within-family id columns
        colIDs = famFile.loc[(famFile[1] == donorID)][[0, 1]]
        removePLINK.append(colIDs)
    else:
        print(f"Warning: No match found for {row.Donor} in {row.Batch}")

if removePLINK:
    removePLINK_df = pd.concat(removePLINK)
    removePLINK_df.to_csv("remove.plink", sep=" ", header=False, index=False)
else:
    # Create empty file if no samples to remove
    with open("remove.plink", "w") as f:
        pass
