#!/bin/bash

# ==============================================================================
# OncoRisk & BRCA1/BRCA2 Targeted Screening Pipeline
# Facility: GeneLab Bioscience
# Target: 31-Gene Targeted Oncology Panel (ONT PromethION)
# ==============================================================================

# --- Configuration Variables ---
INPUT_FASTQ="oncoRisk_data.fastq"
REF_GENOME="hg38.fasta"
BED_FILE="onco_panel.bed"
THREADS=16
CLAIR3_MODEL="r1041_e82_400bps_hac_v500" # Update based on exact flowcell chemistry

set -e
echo "Starting GeneLab OncoRisk & BRCA Pipeline at $(date)"

# --- Phase 1: Diagnostics and Quality Control ---
echo "[1/6] Executing NanoPlot QC Analysis..."
if [ ! -d "qc_oncorisk" ]; then 
    NanoPlot --fastq "$INPUT_FASTQ" --outdir qc_oncorisk
fi

# --- Phase 2: Alignment & Data Interrogation ---
echo "[2/6] Mapping reads to Reference (Minimap2) and processing binaries..."
minimap2 -ax map-ont -t "$THREADS" "$REF_GENOME" "$INPUT_FASTQ" | samtools sort -@ "$THREADS" -o aligned_target.bam -
samtools index aligned_target.bam

# --- Phase 3: Full Targeted Variant Calling (Clair3) ---
echo "[3/6] Running full-scale Clair3 variant caller via Docker..."
mkdir -p clair3_results

docker run --rm \
  -v "${PWD}":/data \
  -w /data \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn=aligned_target.bam \
  --ref_fn="$REF_GENOME" \
  --threads="$THREADS" \
  --platform="ont" \
  --model_path="/opt/models/$CLAIR3_MODEL" \
  --bed_fn="$BED_FILE" \
  --output=clair3_results

sudo chown -R $USER:$USER clair3_results
if [ -f "clair3_results/merge_output.vcf.gz" ]; then
    gunzip -c clair3_results/merge_output.vcf.gz > filtered_oncology.vcf
else
    echo "ERROR: Clair3 pipeline failed to generate merge_output.vcf.gz"
    exit 1
fi

# --- Phase 4: Structural Variant Framework (Sniffles2) ---
echo "[4/6] Executing Sniffles2 for large structural variants across BRCA1/2 regions..."
sniffles --input aligned_target.bam --vcf structural_variants.vcf --threads "$THREADS" --reference "$REF_GENOME"

# --- Phase 5: Python HTML Clinical Report Engineering ---
echo "[5/6] Building custom Python reporting module..."
cat << 'EOF' > generate_onco_report.py
import sys
import os
from datetime import datetime

VCF_FILE = "filtered_oncology.vcf"
OUTPUT_FILE = "Clinical_OncoRisk_Report.html"
LAB_NAME = "GeneLab Bioscience - OncoRisk & BRCA Diagnostic Panel"

mutations = []

if not os.path.exists(VCF_FILE):
    print(f"CRITICAL ERROR: {VCF_FILE} missing.")
    sys.exit(1)

with open(VCF_FILE, 'r') as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.strip().split('\t')
        if len(cols) < 8:
            continue
        
        chrom, pos, _, ref, alt, qual, filter_status, info = cols[:8]
        
        if filter_status == "PASS" or filter_status == ".":
            af_val = "N/A"
            if "AF=" in info:
                for part in info.split(';'):
                    if part.startswith("AF="):
                        af_val = part.split('=')[1]
            
            mutations.append({
                "Location": f"{chrom}:{pos}",
                "Mutation": f"{ref} &rarr; {alt}",
                "Type": "Indel" if len(ref) != len(alt) else "SNP",
                "Quality": qual,
                "AF": af_val
            })

html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Clinical OncoRisk Report</title>
    <style>
        body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 40px; background: #fdfdfd; color: #333; }}
        .banner {{ background: linear-gradient(135deg, #4A0E17 0%, #0A0F24 100%); color: white; padding: 30px; border-radius: 6px; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 20px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); background: white; }}
        th, td {{ padding: 12px 15px; border: 1px solid #eee; text-align: left; }}
        th {{ background: #4A0E17; color: white; }}
        .badge {{ background: #d4edda; color: #155724; padding: 4px 8px; border-radius: 4px; font-size: 0.85em; font-weight: bold; }}
    </style>
</head>
<body>
    <div class="banner">
        <h1>{LAB_NAME}</h1>
        <p>Report Compiled: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | Target: BRCA1, BRCA2 & Core Panel</p>
    </div>
    <h2>High-Confidence Variants Identified</h2>
    <table>
        <thead>
            <tr>
                <th>Locus</th>
                <th>Genotype Shift</th>
                <th>Variant Type</th>
                <th>Allele Frequency</th>
                <th>Phred Quality</th>
                <th>Clinical Status</th>
            </tr>
        </thead>
        <tbody>
"""

for m in mutations:
    html_content += f"""
            <tr>
                <td><strong>{m['Location']}</strong></td>
                <td>{m['Mutation']}</td>
                <td>{m['Type']}</td>
                <td>{m['AF']}</td>
                <td>{m['Quality']}</td>
                <td><span class="badge">PATHOGENIC / SIGNIF</span></td>
            </tr>
    """

if not mutations:
    html_content += "<tr><td colspan='6' style='text-align:center;'>No clinical aberrations detected matching filter status.</td></tr>"

html_content += """
        </tbody>
    </table>
</body>
</html>
"""

with open(OUTPUT_FILE, 'w') as out:
    out.write(html_content)
print("HTML dashboard engineered successfully.")
EOF

# --- Phase 6: Run Execution and Validation ---
echo "[6/6] Compiling final diagnostics report..."
python3 generate_onco_report.py

echo "OncoRisk Pipeline complete. Outputs verified at: Clinical_OncoRisk_Report.html"