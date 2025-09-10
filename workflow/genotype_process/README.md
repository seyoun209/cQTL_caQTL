# caQTL Genotype Processing Pipeline

A Snakemake pipeline for quality control and preparation of genotyping data for imputation, ancestry analysis with EIGENSTRAT, and post-imputation QC of genotype VCF files.

## Overview

This pipeline processes raw genotyping data through the following steps:
1. Convert PLINK files to binary format
2. Sample and variant quality control
3. Strand flipping and chromosome renaming for hg38 compatibility
4. VCF preparation for imputation
5. Ancestry analysis using 1000 Genomes reference panel
6. Post-imputation VCF processing and QC


### Software Dependencies
- Python 3.12.4
- PLINK 1.90b6.21
- R 4.5.0
- BCFtools 1.22
- SAMtools 1.22
- VCFtools
- Eigensoft 8.0.0
- snpflip (optional, for strand flipping)

### Reference Data
- 1000 Genomes Phase 3 reference panel (GRCh38)
- Human genome reference sequence (GRCh38)


### 1. Clone Repository
```bash
git clone git@github.com:seyoun209/cQTL_caQTL.git
cd cQTL_caQTL
```

### 2. Install Eigensoft (Optional)
```bash
mkdir -p tools/dist
cd tools/dist
wget https://github.com/DReichLab/EIG/archive/refs/tags/v8.0.0.tar.gz
tar -xf v8.0.0.tar.gz
cd EIG-8.0.0/src
# Check for OpenBLAS
rpm -qa | grep openblas

# Compile
make clobber
make -j 24
make install
cd ..
mv bin ../../.
cd ../../..
```

### 3. Install snpflip (Optional)
```bash
#!/bin/bash
# Script to install snpflip with Python 3.12 compatibility fixes

echo "Installing snpflip with Python 3.12 compatibility fixes..."

# Create and enter tools directory
mkdir -p tools
cd tools

# Clone snpflip if not already present
if [ ! -d "snpflip" ]; then
    echo "Cloning snpflip repository..."
    git clone https://github.com/biocore-ntnu/snpflip.git
fi

cd snpflip

# Fix the license issue in setup.py
echo "Fixing setup.py license compatibility issue..."
sed -i "s/license=\[.*\]/license='GPL-3.0'/" setup.py

# Install snpflip without dependencies
echo "Installing snpflip in development mode..."
pip install -e . --no-deps

# Install compatible dependencies
echo "Installing compatible dependencies..."
pip install docopt natsort pyfaidx

# Clear bash command cache
hash -r

# Test installation
echo "Testing snpflip installation..."
if command -v snpflip &> /dev/null; then
    echo "✓ snpflip successfully installed"
    echo "Location: $(which snpflip)"
else
    echo "✗ snpflip installation failed"
    echo "Note: Pipeline will skip strand flipping if snpflip unavailable"
fi

echo "Installation complete!"
cd ../../..
```

### 4. Set up Virtual Environment
```bash
python3 -m venv env
source env/bin/activate
pip install -r config/requirements.txt
```

### 1. Sample Sheet (`geno.csv`)
Create a comma-separated file with genotyping batch information:
```csv
Proj,Batch,File_Prefix,Genotyping_Directory
CQTL,COA9_10,COA9_10,/path/to/genotyping/data
```

### 2. Configuration File (`config/Geno_config.yaml`)
Edit parameters for your analysis:
```yaml
# Sample information
geno: 'geno.csv'
remove: 'remove.txt'  # List of samples to remove (can be empty)

# QC thresholds
miss_call: 0.1    # Missing rate per SNP (10% missing = 90% call rate)
maf: 0.01         # Minor allele frequency threshold
hwe: .000001      # Hardy-Weinberg equilibrium p-value

# Reference data paths (GRCh38)
ref: '/proj/phanstiel_lab/References/genomes/1000G/GRCh38/ALL/all_snps/1000G.GRCh38.20190312_chrAll_id'
panel: '/proj/phanstiel_lab/References/genomes/1000G/GRCh38/1000G_phase3.Updatedpanel'
sequence: '/users/s/e/seyoun/Ref/genome/GRCh38.primary_assembly.genome.fa'

# Output directory
paths:
  genotype_outdir: geno_output

# Software versions
python: "3.12.4"
plink: "1.90b6.21"
samtools: "1.22"
r: "4.5.0"
eigensoft: "8.0.0"
```

### 3. Chromosome Mapping File (`addchr.txt`)
Create a file to add "chr" prefix for hg38 compatibility:
```
1	chr1
2	chr2
3	chr3
...
22	chr22
X	chrX
Y	chrY
```

### Pre-imputation Processing
```bash
# Submit job to SLURM
sbatch 01_genoPipe_Run.sbatch

# Or run directly with Snakemake
snakemake -s workflow/genotype_process/00_geno_filter.smk \
  --configfile config/Geno_config.yaml \
  --profile config/profile_slurm_geno \
  -j 120 --latency-wait 500
```

### Key Output Files
- `geno_output/imputation/CQTL_COA9_10_chr*_hg38.vcf.gz`: VCF files ready for imputation
- `geno_output/ancestry/CQTL_COA9_10_ancestry.pdf`: Ancestry PCA plot
- `geno_output/reports/CQTL_COA9_10.sexcheck`: Sex chromosome check results

## Imputation with TOPMed Server

### Upload Settings
Use the [TOPMed Imputation Server](https://imputation.biodatacatalyst.nhlbi.nih.gov/#!run/imputationserver%402.0.0-beta3) with these recommended settings:

**Build:** hg38 (matches your output)
**Reference Panel:** TOPMed r3  
**rsq Filter:** 0.3 (excludes poorly imputed variants)  
**Phasing:** Eagle v2.4 (for unphased data)  
**Population:** 
- If your 21 samples are from mixed populations → Select **"Skip"**
- If from a specific population → Select appropriate population vs TOPMed panel
**Mode:** Quality Control & Imputation  
**Meta-imputation file:** Generally not needed for standard analyses

### Post-imputation Processing
After downloading imputed results:

1. Create VCF list file (`vcf_list_hg38.csv`):
```csv
chr,vcf_path
1,/path/to/chr1.dose.vcf.gz
2,/path/to/chr2.dose.vcf.gz
...
22,/path/to/chr22.dose.vcf.gz
```

2. Run post-imputation workflow:
```bash
sbatch 02_genoPipe_afterImputation_Run.sbatch
```

## Troubleshooting

### Common Issues

**PLINK version conflicts:**
- Use PLINK 1.9 (`plink/1.90b6.21`), not PLINK 2.0
- Check available versions: `module avail plink`

**Chromosome naming issues:**
- Ensure `addchr.txt` uses tab separation
- Check that BCFtools successfully adds "chr" prefix

**snpflip installation problems:**
- Pipeline works without snpflip (skips strand flipping)
- Set `snpflip_available: false` in config if installation fails

### Unlocking Failed Jobs
```bash
source env/bin/activate
snakemake -s workflow/genotype_process/00_geno_filter.smk \
  --configfile config/Geno_config.yaml \
  --unlock
```

## File Structure
```
cQTL_caQTL/
├── config/
│   ├── Geno_config.yaml
│   └── profile_slurm_geno/
├── geno_output/
│   ├── binary/
│   ├── merge/
│   ├── ancestry/
│   ├── imputation/
│   └── reports/
├── scripts/geno_workflow/
├── tools/
└── workflow/genotype_process/
```


