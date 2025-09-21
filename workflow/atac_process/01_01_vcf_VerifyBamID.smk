#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VCF Processing Pipeline - CLEAN VERSION
"""

import pandas as pd
import numpy as np
import os, re, subprocess
from pathlib import Path

#------------------------------------------------------
# Helper
def OUT(*parts):
    return os.path.join(config["paths"]["outdir"], *[p for p in parts if p])

def get_vcf_for_sample(wc):
    sample = wc.sampleName
    for group_name, mask in GROUPS.items():
        if sample in samples_df.loc[mask, 'full_sample_name'].values:
            return OUT("geno", group_name, f"{group_name}.vcf.gz")
    raise ValueError(f"Sample {sample} not found in any group")

#------------------------------------------------------
# Load sample sheet
bam_samplesheet_path = os.path.join(config["paths"]["outdir"], config["bam_samplesheet"])
bam_samples = pd.read_csv(bam_samplesheet_path, sep='\t')
samples_df = bam_samples.copy()
s_config = config['sample_logic']
donor_col = s_config['femur_donor_col']

#------------------------------------------------------
# Handle replicates & femur naming
if str(s_config.get('Sample_handing', 'false')).lower() == 'true':
    print("INFO: Sample_handing is 'true'. Applying special replicate/femur logic.")

    rep_col = s_config['replicate_column']
    is_replicate = samples_df[rep_col].str.contains('replicate', case=False, na=False)

    tissue_col = s_config['femur_column']
    is_femur = samples_df[tissue_col].str.contains('Femur', na=False) | \
               samples_df[donor_col].str.endswith('FE', na=False)

    # Full name (used in BAMs)
    full_df = samples_df.copy()
    full_df['replicate_suffix'] = np.where(is_replicate, 'replicate', '1')
    full_df[donor_col] = np.where(
        is_femur & ~full_df[donor_col].str.endswith('FE'),
        full_df[donor_col] + 'FE',
        full_df[donor_col]
    )
    full_name_cols = s_config['base_name_columns'] + ['replicate_suffix']
    samples_df['full_sample_name'] = full_df[full_name_cols].astype(str).agg('_'.join, axis=1)

    # Standard name (used for main donors only)
    std_df = samples_df.copy()
    std_df['suffix'] = '1'
    std_name_cols = s_config['base_name_columns'] + ['suffix']
    samples_df['standard_sample_name'] = std_df[std_name_cols].astype(str).agg('_'.join, axis=1)

    is_special = is_replicate | is_femur
else:
    print("INFO: Sample_handing is 'false'. Using standard naming only.")
    std_name = samples_df[s_config['base_name_columns']].astype(str).agg('_'.join, axis=1) + '_1'
    samples_df['full_sample_name'] = std_name
    samples_df['standard_sample_name'] = std_name
    is_special = pd.Series(False, index=samples_df.index)

#------------------------------------------------------
# Group donors into 4 categories
group_col = s_config['condition_column']
is_pbs = samples_df[group_col] == s_config['control_group_value']
is_fnf = samples_df[group_col] == s_config['treatment_group_value']

GROUPS = {
    "pbs_geno": is_pbs & ~is_special,
    "fnf_geno": is_fnf & ~is_special,
    "pbs_special": is_pbs & is_special,
    "fnf_special": is_fnf & is_special
}
GROUPS = {name: mask for name, mask in GROUPS.items() if mask.any()}

ALL_SAMPLES = samples_df['full_sample_name'].tolist()

print(f"INFO: {len(ALL_SAMPLES)} total samples -> {len(GROUPS)} groups.")

# Getting the bam
BAM_FILES = {}
BAI_FILES = {}
for _, row in samples_df.iterrows():
    if pd.notna(row.get('Bam_directory')) and pd.notna(row.get('Bam_file')):
        sample_name = row['full_sample_name']
        BAM_FILES[sample_name] = os.path.join(row['Bam_directory'], row['Bam_file'])
        BAI_FILES[sample_name] = BAM_FILES[sample_name] + ".bai"

print(f"INFO: {len(BAM_FILES)} samples have BAM files for verifyBamID")
#------------------------------------------------------
rule all:
    input:
        expand(OUT("geno", "sample_lists", "{group}_samples.txt"), group=GROUPS.keys()),
        expand(OUT("geno", "sample_lists", "{group}_mapping.txt"), group=GROUPS.keys()),
        expand(OUT("geno", "{group}", "{group}.vcf.gz"), group=GROUPS.keys()),
        expand(OUT("QC", "{sampleName}_verifyBamID.selfSM"), sampleName=list(BAM_FILES.keys()))

#------------------------------------------------------
rule write_sample_files:
    input:
        vcf=config["vcf_options"]["original_vcf"]
    output:
        sample_lists=expand(OUT("geno", "sample_lists", "{group}_samples.txt"), group=GROUPS.keys()),
        mapping_files=expand(OUT("geno", "sample_lists", "{group}_mapping.txt"), group=GROUPS.keys())
    params:
        bcftools_version=config["bcftools"]
    run:
        # Get VCF samples
        bcftools_cmd = f"module load bcftools/{params.bcftools_version}; bcftools query -l {input.vcf}"
        result = subprocess.run(bcftools_cmd, shell=True, capture_output=True, text=True, check=True)
        vcf_samples = result.stdout.strip().split('\n')

        output_dir = OUT("geno", "sample_lists")
        os.makedirs(output_dir, exist_ok=True)

        print(f"INFO: Found {len(vcf_samples)} samples in VCF")

        # Create lookup tables
        ctl_lookup = {}
        fnf_lookup = {}
        
        for _, row in samples_df.iterrows():
            donor_id = row[donor_col]
            condition = row[group_col]
            
            if condition == s_config['control_group_value']:
                ctl_lookup[donor_id] = {
                    'standard': row['standard_sample_name'],
                    'full': row['full_sample_name']
                }
            elif condition == s_config['treatment_group_value']:
                fnf_lookup[donor_id] = {
                    'standard': row['standard_sample_name'],
                    'full': row['full_sample_name']
                }

        # Process each group
        for group_name in GROUPS.keys():
            print(f"Processing group: {group_name}")
            
            # Get samples for this group
            group_samples = samples_df.loc[GROUPS[group_name], 'full_sample_name'].tolist()
            
            # Write sample list
            with open(OUT("geno", "sample_lists", f"{group_name}_samples.txt"), "w") as f:
                f.write("\n".join(group_samples))
            
            # Create mapping
            lookup = ctl_lookup if "pbs" in group_name else fnf_lookup
            name_type = "full" if "special" in group_name else "standard"
            
            mappings = []
            remaining_vcf = set(vcf_samples)
            
            # Pass 1: Direct AM matches
            for donor_id, names in lookup.items():
                am_match = re.search(r'AM(\d{4})', donor_id)
                if am_match:
                    am_num = am_match.group(1)
                    for vcf_sample in list(remaining_vcf):
                        if re.match(r'^[0-9]+_AM' + am_num + '$', vcf_sample):
                            mappings.append((vcf_sample, names[name_type]))
                            remaining_vcf.remove(vcf_sample)
                            print(f"Direct match: {vcf_sample} -> {names[name_type]}")
                            break
            
            # Pass 2: Phanstiel matches
            if "special" in group_name:
                # Special groups: try femur samples first, then ankle
                femur_donors = [donor_id for donor_id in lookup.keys() if donor_id.endswith('FE')]
                ankle_donors = [donor_id for donor_id in lookup.keys() if not donor_id.endswith('FE')]
                donor_priority = [femur_donors, ankle_donors]
            else:
                # Regular groups: try ankle samples first, then femur
                ankle_donors = [donor_id for donor_id in lookup.keys() if not donor_id.endswith('FE')]
                femur_donors = [donor_id for donor_id in lookup.keys() if donor_id.endswith('FE')]
                donor_priority = [ankle_donors, femur_donors]
            
            for donor_list in donor_priority:
                for donor_id in donor_list:
                    names = lookup[donor_id]
                    # Skip if already mapped
                    if names[name_type] in [m[1] for m in mappings]:
                        continue
                        
                    am_match = re.search(r'AM(\d{4})', donor_id)
                    if am_match:
                        base_am_num = am_match.group(1)
                        last_two = base_am_num[-2:]
                        
                        for vcf_sample in list(remaining_vcf):
                            if re.match(r'^[0-9]+_Phanstiel_' + last_two + '$', vcf_sample):
                                mappings.append((vcf_sample, names[name_type]))
                                remaining_vcf.remove(vcf_sample)
                                tissue_type = 'femur' if donor_id.endswith('FE') else 'ankle'
                                print(f"Phanstiel match: {vcf_sample} -> {names[name_type]} ({tissue_type})")
                                break
            
            # Add unmapped samples
            for vcf_sample in remaining_vcf:
                mappings.append((vcf_sample, vcf_sample))
            
            # Write mapping file with proper numerical sorting
            def sort_key(mapping_tuple):
                vcf_name = mapping_tuple[0]
                match = re.match(r'^(\d+)_', vcf_name)
                if match:
                    return int(match.group(1))
                return 999
            
            with open(OUT("geno", "sample_lists", f"{group_name}_mapping.txt"), "w") as f:
                for vcf_name, target_name in sorted(mappings, key=sort_key):
                    f.write(f"{vcf_name}\t{target_name}\n")
            
            mapped_count = sum(1 for v, t in mappings if v != t)
            print(f"Group {group_name}: {mapped_count} samples mapped")
#------------------------------------------------------
rule create_subset_vcf:
    input:
        vcf=config["vcf_options"]["original_vcf"],
        mapping_file=OUT("geno", "sample_lists", "{group}_mapping.txt")
    output:
        vcf=OUT("geno", "{group}", "{group}.vcf.gz"),
        tbi=OUT("geno", "{group}", "{group}.vcf.gz.tbi")
    log:
        OUT("logs", "create_subset_vcf_{group}.log")
    params:
        bcftools_version=config["bcftools"],
        samtools_version=config["samtoolsVers"],
    shell:
        """
        module load bcftools/{params.bcftools_version}
        module load samtools/{params.samtools_version}
        mkdir -p $(dirname {output.vcf})

        bcftools reheader --samples {input.mapping_file} {input.vcf} | \
        bcftools view -O z -o {output.vcf} > {log} 2>&1

        tabix -p vcf {output.vcf} >> {log} 2>&1
        """

#------------------------------------------------------
rule verifybamid:
    input:
        bam=lambda wc: BAM_FILES[wc.sampleName],
        bai=lambda wc: BAI_FILES[wc.sampleName],
        vcf=get_vcf_for_sample
    output:
        selfSM=OUT("QC", "{sampleName}_verifyBamID.selfSM"),
        selfRG=OUT("QC", "{sampleName}_verifyBamID.selfRG"),
        bestSM=OUT("QC", "{sampleName}_verifyBamID.bestSM"),
        bestRG=OUT("QC", "{sampleName}_verifyBamID.bestRG"),
        depthSM=OUT("QC", "{sampleName}_verifyBamID.depthSM"),
        depthRG=OUT("QC", "{sampleName}_verifyBamID.depthRG")
    params:
        verifybamid=config['verifybamid'],
        out_prefix=lambda wc: OUT("QC", f"{wc.sampleName}_verifyBamID")
    log:
        OUT("logs", "verifyBamID_{sampleName}.log")
    shell:
        """
        mkdir -p $(dirname {output.selfSM})
        mkdir -p $(dirname {log})
        
        {params.verifybamid} --vcf {input.vcf} --bam {input.bam} --bai {input.bai} \
        --best --out {params.out_prefix} 2> {log}
        """
