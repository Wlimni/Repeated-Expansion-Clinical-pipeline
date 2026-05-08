#!/usr/bin/env bash
# =============================================================================
# Short-read Multi-locus STR Analysis Pipeline with ExpansionHunter + STRipy + REViewer
# ~ 5-10 min runtime (depending on BAM size and number of loci)
# =============================================================================
# Purpose: Detect pathogenic repeat expansions from Illumina short-read WGS data
# Input:   Sample BAM file (e.g., G241713.bam) aligned to hg38
# Output:  
#   - multi_realigned.bam           : Reads realigned to repeat graph
#   - multi_realigned.sorted.bam    : Sorted version for REViewer input
#   - ${SAMPLE}_multi.vcf           : Raw ExpansionHunter VCF output
#   - ${SAMPLE}_multi_annotated.vcf : Annotated with STRipy clinical info
#   - ${SAMPLE}_multi_info.tsv      : Full summary table (all loci)
#   - ${SAMPLE}_significant_only_info.tsv : Filtered table (pathogenic only)
#   - reviewer/*.svg                : REViewer plots for each locus
# =============================================================================
# Dependencies: ExpansionHunter v5.0.0+, REViewer v0.2.7+, samtools, bcftools, curl
# =============================================================================
# conda install -c conda-forge -c bioconda expansionhunter reviewer samtools bcftools -y

set -e
set -u

# =============================================================================
# CONFIGURATION
# =============================================================================

SAMPLE="G241713"
REF="/home/gglab2/wgs/data/ref/Homo_sapiens_assembly38.fasta"
FINAL_BAM="/home/gglab2/wgs/data/align_G241713_SR/G241713.bam"
CATALOG="/home/gglab2/wgs/data/SR/known_variant_catalog.json" # (From STRipy, includes off-target)
OUTDIR="/home/gglab2/wgs/results_SR/${SAMPLE}/EH_multi"

# =============================================================================
# SETUP
# =============================================================================

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "================================================================================"
echo "Short-read STR Analysis Pipeline (ExpansionHunter + STRipy + REViewer)"
echo "Sample: $SAMPLE"
echo "Started: $(date)"
echo "================================================================================"

# =============================================================================
# STEP 1: ExpansionHunter
# =============================================================================
echo ""
echo "[Step 1/6] Running ExpansionHunter on multi-locus catalog"

if [[ -s "${SAMPLE}_multi.vcf" ]]; then
    echo "  ✓ VCF already exists → skipping ExpansionHunter"
else
    echo "  → Running ExpansionHunter..."
    ExpansionHunter \
        --reads "$FINAL_BAM" \
        --reference "$REF" \
        --variant-catalog "$CATALOG" \
        --output-prefix "${SAMPLE}_multi" \
        --threads 12 \
        --analysis-mode streaming
    echo "  ✓ ExpansionHunter completed"
fi

# =============================================================================
# STEP 2: Prepare BAMlet for REViewer
# =============================================================================
echo ""
echo "[Step 2/6] Preparing BAMlet for REViewer"

if [[ -f "${SAMPLE}_multi_realigned.sorted.bam" ]]; then
    echo "  ✓ Sorted BAMlet already exists"
else
    echo "  → Sorting and indexing BAMlet..."
    samtools sort -@ 8 -o "${SAMPLE}_multi_realigned.sorted.bam" "${SAMPLE}_multi_realigned.bam"
    samtools index "${SAMPLE}_multi_realigned.sorted.bam"
    echo "  ✓ BAMlet ready"
fi

# =============================================================================
# STEP 3: STRipy Annotation
# =============================================================================
echo ""
echo "[Step 3/6] Annotating VCF with STRipy"

ANNOTATED_VCF="${SAMPLE}_multi_annotated.vcf"

if [[ -s "$ANNOTATED_VCF" ]]; then
    echo "  ✓ Annotated VCF already exists"
else
    echo "  → Sending VCF to STRipy API..."
    curl -F "file=@${SAMPLE}_multi.vcf" https://api.stripy.org/annotateVCF > "$ANNOTATED_VCF"
    echo "  ✓ STRipy annotation completed"
fi

# =============================================================================
# STEP 4: Create Summary Table (Robust VCF parsing)
# =============================================================================
echo ""
echo "[Step 4/6] Creating Summary Table"

SUMMARY_FILE="${SAMPLE}_multi_info.tsv"

# Create header with proper column names
cat > "$SUMMARY_FILE" << 'EOF'
#CHROM	POS	FILTER	END	Reference_Copy_Number	Reference_Length_bp	Repeat_Unit	Variant_ID	Disease_Name	Inheritance	Clinical_Range	Genotype	Support_Type	Repeat_Count	Confidence_Interval	Spanning_Reads	Flanking_Reads	InRepeat_Reads	Locus_Coverage
EOF

# Parse VCF with robust extraction
awk -F'\t' '!/^#/ && NF >= 10 {
    # Skip if CHROM or POS is empty
    if ($1 == "" || $2 == "") next;
    
    # ----- Extract INFO fields (robust method) -----
    varid = ""; disname = ""; disinher = ""; disrange = "";
    end = ""; ref_cn = ""; rl = ""; ru = "";
    
    split($8, info, ";");
    for(i in info) {
        # Find position of "=" and extract after it
        eq_pos = index(info[i], "=");
        if(eq_pos > 0) {
            key = substr(info[i], 1, eq_pos - 1);
            val = substr(info[i], eq_pos + 1);
            
            if(key == "END") end = val;
            else if(key == "REF") ref_cn = val;
            else if(key == "RL") rl = val;
            else if(key == "RU") ru = val;
            else if(key == "VARID") varid = val;
            else if(key == "DISNAME") disname = val;
            else if(key == "DISINHER") disinher = val;
            else if(key == "DISRANGE") disrange = val;
        }
    }
    
    # Skip if no Variant_ID (likely not a valid repeat locus)
    if (varid == "") next;
    
    # ----- Extract FORMAT sample data -----
    split($9, fmt, ":");
    split($10, samp, ":");
    
    # Default values
    gt = "./."; so = "./."; repcn = "./."; repci = "./."; 
    adsp = "./."; adfl = "./."; adir = "./."; lc = "0";
    
    # Map format to values
    for(i in fmt) {
        val = samp[i];
        if(fmt[i] == "GT") gt = val;
        else if(fmt[i] == "SO") so = val;
        else if(fmt[i] == "REPCN") repcn = val;
        else if(fmt[i] == "REPCI") repci = val;
        else if(fmt[i] == "ADSP") adsp = val;
        else if(fmt[i] == "ADFL") adfl = val;
        else if(fmt[i] == "ADIR") adir = val;
        else if(fmt[i] == "LC") lc = val;
    }
    
    # Skip rows where Genotype is "./." (no call) AND Repeat_Count is "./."
    if (gt == "./." && repcn == "./.") next;
    
    # ----- Output -----
    print $1 "\t" $2 "\t" $7 "\t" end "\t" ref_cn "\t" rl "\t" ru "\t" varid "\t" disname "\t" disinher "\t" disrange "\t" gt "\t" so "\t" repcn "\t" repci "\t" adsp "\t" adfl "\t" adir "\t" lc;
}' "$ANNOTATED_VCF" >> "$SUMMARY_FILE"

# Remove any empty lines from the file
sed -i '/^$/d' "$SUMMARY_FILE"

# Verify and preview
if [[ -s "$SUMMARY_FILE" ]]; then
    # Get actual row count (excluding header)
    row_count=$(tail -n +2 "$SUMMARY_FILE" | wc -l)
    
    echo "  ✓ Summary table created: $SUMMARY_FILE"
    echo "  ✓ Total rows (excluding header): $row_count"
    echo ""
    echo "  Preview (first 5 data rows):"
    echo "  ---------------------------------------------------"
    head -6 "$SUMMARY_FILE" | cut -f1,2,4,5,6,7,8,9,10,11
    echo "  ---------------------------------------------------"
    
    # Quick validation check
    echo ""
    echo "  Validation check (Inheritance column):"
    tail -n +2 "$SUMMARY_FILE" | cut -f10 | sort | uniq -c
else
    echo "  ❌ ERROR: Failed to create summary table"
    exit 1
fi

# =============================================================================
# STEP 5: Final Summary
# =============================================================================
echo ""
echo "[Step 5/6] Final Summary"

# Count total loci (skip header lines)
total=$(grep -v "^#" "$ANNOTATED_VCF" | grep -c -v "^$" || echo "0")

# Count significant loci (Pathogenic or Intermediate in DISRANGE)
significant=$(grep -v "^#" "$ANNOTATED_VCF" | grep -c -E "DISRANGE=.*(Pathogenic|Intermediate)" || echo "0")

echo "  ┌─────────────────────────────────────────┐"
echo "  │  SAMPLE: $SAMPLE"
echo "  ├─────────────────────────────────────────┤"
echo "  │  Total loci analyzed:       $total"
echo "  │  Potentially significant:   $significant"
echo "  └─────────────────────────────────────────┘"

if [[ $significant -gt 0 ]]; then
    echo ""
    echo "  ⚠️  CLINICAL ALERT: Potentially significant loci found:"
    grep -v "^#" "$ANNOTATED_VCF" | grep -E "DISRANGE=.*(Pathogenic|Intermediate)" | while read -r line; do
        varid=$(echo "$line" | grep -oP 'VARID=\K[^;]+' || echo "UNKNOWN")
        disrange=$(echo "$line" | grep -oP 'DISRANGE=\K[^;]+' || echo "N/A")
        echo "    • $varid: $disrange"
    done
fi

# =============================================================================
# STEP 5.5: Create significant-only report
# =============================================================================
echo ""
echo "[Step 5.5/6] Creating significant-only report"

SIG_REPORT="${SAMPLE}_significant_only_info.tsv"

if [[ $significant -gt 0 ]]; then
    # Create header with same column names
    cat > "$SIG_REPORT" << 'EOF'
#CHROM	POS	FILTER	END	Reference_Copy_Number	Reference_Length_bp	Repeat_Unit	Variant_ID	Disease_Name	Inheritance	Clinical_Range	Genotype	Support_Type	Repeat_Count	Confidence_Interval	Spanning_Reads	Flanking_Reads	InRepeat_Reads	Locus_Coverage
EOF
    
    # Parse only significant lines (Pathogenic or Intermediate in DISRANGE)
    awk -F'\t' '!/^#/ && $8 ~ /Pathogenic|Intermediate/ {
        # ----- Extract INFO fields -----
        varid = ""; disname = ""; disinher = ""; disrange = "";
        end = ""; ref_cn = ""; rl = ""; ru = "";
        
        split($8, info, ";");
        for(i in info) {
            eq_pos = index(info[i], "=");
            if(eq_pos > 0) {
                key = substr(info[i], 1, eq_pos - 1);
                val = substr(info[i], eq_pos + 1);
                
                if(key == "END") end = val;
                else if(key == "REF") ref_cn = val;
                else if(key == "RL") rl = val;
                else if(key == "RU") ru = val;
                else if(key == "VARID") varid = val;
                else if(key == "DISNAME") disname = val;
                else if(key == "DISINHER") disinher = val;
                else if(key == "DISRANGE") disrange = val;
            }
        }
        
        # ----- Extract FORMAT sample data -----
        split($9, fmt, ":");
        split($10, samp, ":");
        
        gt = "./."; so = "./."; repcn = "./."; repci = "./."; 
        adsp = "./."; adfl = "./."; adir = "./."; lc = "0";
        
        for(i in fmt) {
            val = samp[i];
            if(fmt[i] == "GT") gt = val;
            else if(fmt[i] == "SO") so = val;
            else if(fmt[i] == "REPCN") repcn = val;
            else if(fmt[i] == "REPCI") repci = val;
            else if(fmt[i] == "ADSP") adsp = val;
            else if(fmt[i] == "ADFL") adfl = val;
            else if(fmt[i] == "ADIR") adir = val;
            else if(fmt[i] == "LC") lc = val;
        }
        
        # ----- Output -----
        print $1 "\t" $2 "\t" $7 "\t" end "\t" ref_cn "\t" rl "\t" ru "\t" varid "\t" disname "\t" disinher "\t" disrange "\t" gt "\t" so "\t" repcn "\t" repci "\t" adsp "\t" adfl "\t" adir "\t" lc;
    }' "$ANNOTATED_VCF" >> "$SIG_REPORT"
    
    # Remove empty lines
    sed -i '/^$/d' "$SIG_REPORT"
    
    sig_row_count=$(tail -n +2 "$SIG_REPORT" | wc -l)
    echo "  ✓ Significant-only report: $SIG_REPORT ($sig_row_count loci)"
    echo ""
    echo "  Significant loci summary:"
    echo "  -------------------------"
    tail -n +2 "$SIG_REPORT" | cut -f8,10,11,14 | head -10
    echo "  -------------------------"
else
    echo "  ✅ No significant loci found - skipping significant-only report"
fi

# =============================================================================
# STEP 6: Plot significant loci with REViewer
# =============================================================================
echo ""
echo "[Step 6/6] Plotting significant loci with REViewer"

PLOTS_DIR="${OUTDIR}/reviewer"
mkdir -p "$PLOTS_DIR"

if [[ $significant -gt 0 ]]; then
    # Extract unique locus names from significant lines
    grep -v "^#" "$ANNOTATED_VCF" | grep -E "DISRANGE=.*(Pathogenic|Intermediate)" | grep -oP 'VARID=\K[^;]+' | sort -u | while read -r LOCUS; do
        echo "  → Plotting: $LOCUS"
        
        REViewer \
            --reads "${SAMPLE}_multi_realigned.sorted.bam" \
            --vcf "$ANNOTATED_VCF" \
            --reference "$REF" \
            --catalog "$CATALOG" \
            --locus "$LOCUS" \
            --output-prefix "${PLOTS_DIR}/${LOCUS}" || true
        
        # Rename plot files to consistent naming
        if [[ -f "${PLOTS_DIR}/${LOCUS}.${LOCUS}.svg" ]]; then
            mv "${PLOTS_DIR}/${LOCUS}.${LOCUS}.svg" "${PLOTS_DIR}/${SAMPLE}_reviewer_${LOCUS}.svg"
            echo "    ✓ Plot saved: ${PLOTS_DIR}/${SAMPLE}_reviewer_${LOCUS}.svg"
        elif [[ -f "${PLOTS_DIR}/${LOCUS}.svg" ]]; then
            mv "${PLOTS_DIR}/${LOCUS}.svg" "${PLOTS_DIR}/${SAMPLE}_reviewer_${LOCUS}.svg"
            echo "    ✓ Plot saved: ${PLOTS_DIR}/${SAMPLE}_reviewer_${LOCUS}.svg"
        else
            echo "    ⚠️  Plot not found for ${LOCUS}"
        fi
    done
else
    echo "  ✅ No significant loci - skipping REViewer plots"
fi

echo "  ✓ Plotting complete"

# =============================================================================
# CLEANUP: Delete all empty log files
# =============================================================================
echo ""
echo "[Cleanup] Removing empty log files..."

find "$OUTDIR" -type f -name "*.log" -size 0 -delete 2>/dev/null

echo "  ✓ Empty log files removed"

# =============================================================================
# PIPELINE COMPLETE
# =============================================================================
echo ""
echo "================================================================================"
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "================================================================================"
echo "Sample: $SAMPLE"
echo "Summary table: $SUMMARY_FILE"
echo "Significant-only: $SIG_REPORT"
echo "Plots directory: $PLOTS_DIR"
echo "Completion time: $(date)"
echo "================================================================================"

# List output files
echo ""
echo "Output files:"
ls -lh "$SUMMARY_FILE" 2>/dev/null || echo "  No summary file"
ls -lh "$SIG_REPORT" 2>/dev/null || echo "  No significant report"
ls -lh "$PLOTS_DIR"/*.svg 2>/dev/null || echo "  No plot files found"
