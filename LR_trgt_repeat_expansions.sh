#!/usr/bin/env bash
# =============================================================================
# Long-read Multi-locus STR Analysis Pipeline with TRGT + TRVZ Analysis
# ~ 10s in this example
# =============================================================================
# Purpose: Detect pathogenic repeat expansions from PacBio long-read data
# Input:  Sample IDs (SAMPLES)
# Output: Genotype summary table + significant-only report + plots for significant repeats
# =============================================================================

set -e  # Stop on error
set -u  # Error on undefined variables

# =============================================================================
# USER CONFIGURATION - EDIT THESE AS NEEDED
# =============================================================================

# List of sample IDs to analyze
SAMPLES=("242979" "241819" "241713")

# Reference genome (GRCh38/hg38)
REFERENCE="/home/gglab2/wgs/data/ref/Homo_sapiens_assembly38.fasta"

# Bed file with pathogenic repeat regions (56 loci) (from https://github.com/PacificBiosciences/trgt/blob/main/repeats/pathogenic_repeats.hg38.bed)
REPEAT_CATALOG="/home/gglab2/wgs/data/pathogenic_repeats.hg38.bed"

# Main output directory
OUTPUT_BASE="/home/gglab2/wgs/results_LR"

# Premutation thresholds per gene (cols: gene, premutation_threshold, pathogenic_motif) (from https://github.com/PacificBiosciences/TRexs/blob/main/resources/repeats_information.tsv)
THRESHOLD_FILE="/home/gglab2/wgs/data/premutation_thresholds.tsv"

# Path to PacBio BAM files (pattern: /path/${SAMPLE}/alignment/${SAMPLE}.GRCh38.haplotagged.bam)
BAM_BASE="/home/gglab2/share/lrWGS/PacBio"

# =============================================================================
# SETUP
# =============================================================================

cd /home/gglab2 || { echo "ERROR: Cannot access /home/gglab2"; exit 1; }
mkdir -p "$OUTPUT_BASE"

# Create timestamp for logging
LOG_FILE="${OUTPUT_BASE}/pipeline_$(date +%Y%m%d_%H%M%S).log"

# Start logging
exec > >(tee -a "$LOG_FILE")
exec 2>&1

echo "================================================================================"
echo "Long-read STR Analysis Pipeline"
echo "Started: $(date)"
echo "================================================================================"
echo "Samples to process: ${SAMPLES[@]}"
echo "Output directory: $OUTPUT_BASE"
echo "Log file: $LOG_FILE"
echo "================================================================================"

# =============================================================================
# LOAD THRESHOLDS FOR PATHOGENIC REPEATS
# =============================================================================

[[ -f "$THRESHOLD_FILE" ]] || { echo "ERROR: $THRESHOLD_FILE not found"; exit 1; }

declare -A THRESHOLD
declare -A PATHOGENIC_MOTIF

while IFS=$'\t' read -r gene threshold path_motif; do
    [[ "$gene" == "gene" ]] && continue
    THRESHOLD["$gene"]="$threshold"
    PATHOGENIC_MOTIF["$gene"]="$path_motif"
done < "$THRESHOLD_FILE"

# =============================================================================
# HELPER FUNCTION: Extract pathogenic motif count from TRGT output
# =============================================================================

get_pathogenic_count() {
    local mc_field="$1"           # Format: "count1_count2,count1_count2" (allele1, allele2)
    local motifs_list="$2"        # Comma-separated motifs (e.g., "CAG,CGG,CAG")
    local target_motif="$3"       # Motif to count (e.g., "CAG")
    
    # Helper: Extract largest number from MC field (fallback for unknown motifs)
    get_largest_number() {
        echo "$1" | grep -oE '[0-9]+' | sort -nr | head -1
    }
    
    # Case 1: No specific pathogenic motif defined → use largest count
    if [[ -z "$target_motif" || "$target_motif" == "N/A" || "$target_motif" == "." ]]; then
        get_largest_number "$mc_field"
        return
    fi
    
    # Case 2: Find which position the target motif occupies in the MOTIFS list
    local motif_index=-1
    IFS=',' read -ra motifs <<< "$motifs_list"
    for i in "${!motifs[@]}"; do
        if [[ "${motifs[$i]}" == "$target_motif" ]]; then
            motif_index=$i
            break
        fi
    done
    
    # Case 3: Motif not found → fallback to largest number
    if [[ $motif_index -eq -1 ]]; then
        get_largest_number "$mc_field"
        return
    fi
    
    # Case 4: Extract count for specific motif from both alleles
    local position=$((motif_index + 1))
    local allele1=$(echo "$mc_field" | cut -d',' -f1 | cut -d'_' -f$position | grep -o '^[0-9]\+' || echo "0")
    local allele2=$(echo "$mc_field" | cut -d',' -f2 | cut -d'_' -f$position | grep -o '^[0-9]\+' || echo "0")
    
    # Return the larger of the two alleles
    if [[ "$allele1" -gt "$allele2" ]]; then
        echo "$allele1"
    else
        echo "$allele2"
    fi
}

# =============================================================================
# MAIN PIPELINE: Process each sample
# =============================================================================

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "================================================================================"
    echo "Processing Sample: $SAMPLE"
    echo "================================================================================"
    
    # Define file paths
    BAM_FILE="${BAM_BASE}/${SAMPLE}/alignment/${SAMPLE}.GRCh38.haplotagged.bam"
    SAMPLE_OUTDIR="${OUTPUT_BASE}/${SAMPLE}/TRGT"
    VCF_FILE="${SAMPLE_OUTDIR}/${SAMPLE}.vcf.gz"
    SPANNING_BAM="${SAMPLE_OUTDIR}/${SAMPLE}.spanning.sorted.bam"
    SUMMARY_FILE="${SAMPLE_OUTDIR}/${SAMPLE}_trgt_summary.tsv"
    
    # Create sample output directory
    mkdir -p "$SAMPLE_OUTDIR"
    cd "$SAMPLE_OUTDIR" || { echo "ERROR: Cannot access $SAMPLE_OUTDIR"; continue; }
    
    # Verify input BAM exists
    if [[ ! -f "$BAM_FILE" ]]; then
        echo "⚠️  SKIPPING $SAMPLE: BAM file not found"
        echo "   Expected: $BAM_FILE"
        continue
    fi
    
    echo "Input BAM: $BAM_FILE"
    echo "Output dir: $SAMPLE_OUTDIR"
    
    # -------------------------------------------------------------------------
    # STEP 1: Run TRGT genotyping
    # -------------------------------------------------------------------------
    echo ""
    echo "[Step 1/5] TRGT Genotyping..."
    
    if [[ -f "$VCF_FILE" ]]; then
        echo "  ✓ VCF already exists: $VCF_FILE"
    else
        echo "  → Running TRGT genotype..."
        trgt genotype \
            --genome "$REFERENCE" \
            --repeats "$REPEAT_CATALOG" \
            --reads "$BAM_FILE" \
            --output-prefix "$SAMPLE" \
            --threads 12 \
            2> "${SAMPLE}_trgt_error.log"
        
        if [[ -f "${SAMPLE}.vcf.gz" ]]; then
            # Sort and index VCF
            bcftools sort -Ob -o "${SAMPLE}.sorted.vcf.gz" "${SAMPLE}.vcf.gz"
            bcftools index "${SAMPLE}.sorted.vcf.gz"
            mv "${SAMPLE}.sorted.vcf.gz" "$VCF_FILE"
            echo "  ✓ TRGT completed successfully"
        else
            echo "  ❌ ERROR: TRGT failed for $SAMPLE"
            continue
        fi
    fi
    
    # -------------------------------------------------------------------------
    # STEP 2: Prepare spanning reads BAM for visualization
    # -------------------------------------------------------------------------
    echo ""
    echo "[Step 2/5] Preparing spanning reads BAM..."
    
    if [[ -f "$SPANNING_BAM" ]]; then
        echo "  ✓ Spanning BAM already exists and sorted"
    elif [[ -f "${SAMPLE}.spanning.bam" ]]; then
        echo "  → Sorting spanning BAM..."
        samtools sort -@ 8 -o "$SPANNING_BAM" "${SAMPLE}.spanning.bam"
        samtools index "$SPANNING_BAM"
        echo "  ✓ Spanning BAM ready"
    else
        echo "  ⚠️  No spanning BAM found (will skip TRVZ plots)"
    fi
    
    # -------------------------------------------------------------------------
    # STEP 3: Create summary table with pathogenicity classification
    # -------------------------------------------------------------------------
    echo ""
    echo "[Step 3/5] Creating summary table..."
    
    {
        # Header
        echo -e "Gene\tGenotype\tMC(MotifCount)\tSD(SpanningDepth)\tAL(AlleleLength)\tAP(Purity)\tAM(Methylation)\tPathogenicCount\tThreshold\tSignificant?"
        
        # Process each variant in VCF
        bcftools view -H "$VCF_FILE" 2>/dev/null | while read -r line; do
            # Extract gene name from TRID
            gene=$(echo "$line" | grep -oP 'TRID=\K[^;]+' || echo "UNKNOWN")
            [[ -z "$gene" || "$gene" == "." ]] && gene="UNKNOWN"
            
            # Extract motifs from INFO field
            motifs=$(echo "$line" | grep -oP 'MOTIFS=\K[^;]+' || echo "")
            
            # Extract genotype and metrics from sample column
            sample_data=$(echo "$line" | cut -f10)
            genotype=$(echo "$sample_data" | cut -d':' -f1)
            allele_length=$(echo "$sample_data" | cut -d':' -f2)
            spanning_depth=$(echo "$sample_data" | cut -d':' -f4)
            motif_counts=$(echo "$sample_data" | cut -d':' -f5)
            purity=$(echo "$sample_data" | cut -d':' -f7)
            methylation=$(echo "$sample_data" | cut -d':' -f8)
            
            # Get threshold and pathogenic motif for this gene
            threshold=${THRESHOLD["$gene"]:-}
            pathogenic_motif=${PATHOGENIC_MOTIF["$gene"]:-""}
            
            # Skip if no threshold defined (gene not in threshold file)
            if [[ -z "$threshold" ]]; then
                continue
            fi
            
            # Calculate pathogenic repeat count
            pathogenic_count=$(get_pathogenic_count "$motif_counts" "$motifs" "$pathogenic_motif")
            
            # Determine significance
            significant="No"
            if [[ "$genotype" != "0/0" ]] && [[ "$pathogenic_count" -ge "$threshold" ]]; then
                significant="Yes"
            fi
            
            # Output line
            echo -e "$gene\t$genotype\t$motif_counts\t$spanning_depth\t$allele_length\t$purity\t$methylation\t$pathogenic_count\t$threshold\t$significant"
        done
    } > "$SUMMARY_FILE"
    
    echo "  ✓ Summary table created: $SUMMARY_FILE"
    echo ""
    
    # -------------------------------------------------------------------------
    # STEP 4: Generate visualizations for significant loci (TRVZ)
    # -------------------------------------------------------------------------
    echo ""
    echo "[Step 4/5] Generating TRVZ plots for significant loci..."
    
    PLOTS_DIR="$SAMPLE_OUTDIR/plots"
    mkdir -p "$PLOTS_DIR"
    
    awk -F'\t' 'NR>1 && $NF=="Yes" {print $1}' "$SUMMARY_FILE" | while read -r gene; do
        gene=$(echo "$gene" | xargs)
        echo "  → Plotting: $gene"
        
        if [[ -f "$SPANNING_BAM" ]]; then
            trgt plot \
                --genome "$REFERENCE" \
                --repeats "$REPEAT_CATALOG" \
                --vcf "$VCF_FILE" \
                --spanning-reads "$SPANNING_BAM" \
                --repeat-id "$gene" \
                --image "${PLOTS_DIR}/${SAMPLE}_${gene}.svg" \
                2> "trvz_${gene}.log"
            
            if [[ -f "${PLOTS_DIR}/${SAMPLE}_${gene}.svg" ]]; then
                echo "    ✓ Plot saved: ${PLOTS_DIR}/${SAMPLE}_${gene}.svg"
            else
                echo "    ⚠️  TRVZ plot failed (see trvz_${gene}.log)"
            fi
        else
            echo "    ⚠️  No spanning BAM available - skipping plot"
        fi
    done
    
    # -------------------------------------------------------------------------
    # STEP 5: Generate clinical summary report
    # -------------------------------------------------------------------------
    echo ""
    echo "[Step 5/5] Generating clinical summary..."
    
    total_loci=$(awk 'NR>1 && $1!="" && $1!~/#/ {print $1}' "$SUMMARY_FILE" | sort -u | wc -l)
    significant_loci=$(awk -F'\t' 'NR>1 && $NF=="Yes" {print $1}' "$SUMMARY_FILE" | sort -u | wc -l)
    
    echo "  ┌─────────────────────────────────────────┐"
    echo "  │  SAMPLE: $SAMPLE"
    echo "  ├─────────────────────────────────────────┤"
    echo "  │  Total loci genotyped:     $total_loci"
    echo "  │  Clinically significant:   $significant_loci"
    echo "  └─────────────────────────────────────────┘"
    
    if [[ $significant_loci -gt 0 ]]; then
        echo ""
        echo "  ⚠️  CLINICAL ALERT: Significant expansions detected:"
        awk -F'\t' 'NR>1 && $NF=="Yes" {printf "    • %s: %s repeats (threshold ≥ %s)\n", $1, $8, $9}' "$SUMMARY_FILE"
        
        echo ""
        echo "  📊 Generating detailed report for significant loci only..."
        
        SIGNIFICANT_REPORT="${SAMPLE_OUTDIR}/${SAMPLE}_significant_only.tsv"
        
        {
            echo -e "Gene\tGenotype\tMC(MotifCount)\tSD(SpanningDepth)\tAL(AlleleLength)\tAP(Purity)\tAM(Methylation)\tPathogenicCount\tThreshold"
            
            awk -F'\t' 'NR>1 && $NF=="Yes" {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' "$SUMMARY_FILE"
        } > "$SIGNIFICANT_REPORT"
        
        echo "  ✓ Significant-only report: $SIGNIFICANT_REPORT"
        echo ""
        echo "  Significant loci details:"
        column -t -s $'\t' "$SIGNIFICANT_REPORT" | head -20
        
    else
        echo "  ✅ No clinically significant expansions found"
    fi
    
    echo ""
    echo "================================================================================"
    echo "Sample $SAMPLE completed successfully"
    echo "================================================================================"
    
done

# =============================================================================
# PIPELINE COMPLETE
# =============================================================================

echo ""
echo "================================================================================"
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "================================================================================"
echo "All samples processed: ${SAMPLES[@]}"
echo "Results saved in: $OUTPUT_BASE"
echo "Log file: $LOG_FILE"
echo "Completion time: $(date)"
echo "================================================================================"
