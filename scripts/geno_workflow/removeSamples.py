
import sys
import pandas as pd
import os

def create_empty_remove_file():
    """Create empty remove.plink file when no samples to remove"""
    with open('remove.plink', 'w') as f:
        pass

def main():
    # Check arguments
    if len(sys.argv) < 2:
        create_empty_remove_file()
        return
    
    remove_file = sys.argv[1]
    
    # Check if remove file exists and has content
    if not os.path.isfile(remove_file):
        create_empty_remove_file()
        return
    
    # Read samples to remove
    try:
        removeSamples = pd.read_csv(remove_file, sep=" ", header=None)
        if removeSamples.empty or removeSamples.shape[1] < 2:
            create_empty_remove_file()
            return
    except:
        create_empty_remove_file()
        return
    
    # Find matching samples in .fam files
    removePLINK = []
    
    # Look for .fam files in binary directories
    for batch_dir in ['geno_output/binary', 'output/binary']:
        if not os.path.exists(batch_dir):
            continue
            
        for file in os.listdir(batch_dir):
            if not file.endswith('.fam'):
                continue
                
            fam_path = os.path.join(batch_dir, file)
            try:
                famFile = pd.read_csv(fam_path, sep=" ", header=None)
                
                # Match samples (assuming remove file has FID IID format)
                for _, remove_row in removeSamples.iterrows():
                    matches = famFile[(famFile[0] == remove_row[0]) & (famFile[1] == remove_row[1])]
                    if not matches.empty:
                        removePLINK.append(matches[[0, 1]])
            except:
                continue
    
    # Write output
    if removePLINK:
        result = pd.concat(removePLINK, ignore_index=True)
        result.to_csv("remove.plink", sep=" ", header=False, index=False)
    else:
        create_empty_remove_file()

if __name__ == "__main__":
    main()
