#!/usr/bin/env bash
# =============================================================================
# Multi-locus STR Analysis with ExpansionHunter + STRipy + REViewer
# ~ 5 min in this example
# =============================================================================
# Purpose: Detect pathogenic repeat expansions from short-read WGS data
# Input:  Single sample ID (edit SAMPLE below)
# Output: Genotype summary table + significant-only report + REViewer plots
# =============================================================================

set -e

# =============================================================================
# CONFIGURATION
# =============================================================================

SAMPLE="G241713"
REF="/home/gglab2/wgs/data/ref/Homo_sapiens_assembly38.fasta"
FINAL_BAM="/home/gglab2/wgs/results_SR/${SAMPLE}/align/${SAMPLE}.bam"
CATALOG="/home/gglab2/wgs/data/known_variant_catalog.json" # (From https://github.com/PacificBiosciences/trgt, include off-target)
OUTDIR="/home/gglab2/wgs/results_SR/${SAMPLE}/EH_multi"

# =============================================================================
# SETUP
# =============================================================================

mkdir -p "$OUTDIR"
cd "$OUTDIR"

EH_PREFIX="${SAMPLE}_multi"

echo "================================================================================"
echo "Multi-locus ExpansionHunter Analysis + STRipy Annotation"
echo "Sample: $SAMPLE"
echo "Started: $(date)"
echo "================================================================================"

# =============================================================================
# STEP 1: ExpansionHunter
# =============================================================================
echo ""
echo "[Step 1/6] Running ExpansionHunter on multi-locus catalog"

if [[ -s "${EH_PREFIX}.vcf" ]]; then
    echo "  ✓ VCF already exists → skipping ExpansionHunter"
else
    echo "  → Running ExpansionHunter..."
    ExpansionHunter \
        --reads "$FINAL_BAM" \
        --reference "$REF" \
        --variant-catalog "$CATALOG" \
        --output-prefix "$EH_PREFIX" \
        --threads 12 \
        --analysis-mode streaming 2> eh.log
    echo "  ✓ ExpansionHunter completed"
fi

# =============================================================================
# STEP 2: Prepare BAMlet
# =============================================================================
echo ""
echo "[Step 2/6] Preparing BAMlet for REViewer"

if [[ -f "${EH_PREFIX}_realigned.sorted.bam" ]]; then
    echo "  ✓ Sorted BAMlet already exists"
else
    echo "  → Sorting and indexing BAMlet..."
    samtools sort -@ 8 -o "${EH_PREFIX}_realigned.sorted.bam" "${EH_PREFIX}_realigned.bam"
    samtools index "${EH_PREFIX}_realigned.sorted.bam"
    echo "  ✓ BAMlet ready"
fi

# =============================================================================
# STEP 3: STRipy Annotation
# =============================================================================
echo ""
echo "[Step 3/6] Annotating VCF with STRipy"

ANNOTATED_VCF="${EH_PREFIX}_annotated.vcf"

if [[ -s "$ANNOTATED_VCF" ]]; then
    echo "  ✓ Annotated VCF already exists"
else
    echo "  → Sending VCF to STRipy API..."
    curl -F "file=@${EH_PREFIX}.vcf" https://api.stripy.org/annotateVCF > "$ANNOTATED_VCF"
    echo "  ✓ STRipy annotation completed"
fi

# =============================================================================
# STEP 4: Create Summary Table
# =============================================================================
echo ""
echo "[Step 4/6] Creating Summary Table"

SUMMARY_FILE="${EH_PREFIX}_summary.tsv"

{
    echo -e "Locus\tGenotype\tRepeatCounts\tDisease\tInheritance\tRange"
    grep -v "^#" "$ANNOTATED_VCF" | while read -r line; do
        # Extract VARID from INFO column (locus name)
        varid=$(echo "$line" | grep -oP 'VARID=\K[^;]+' || echo "UNKNOWN")
        
        # Extract genotype and repeat counts from FORMAT column
        format=$(echo "$line" | cut -f9)
        sample=$(echo "$line" | cut -f10)
        
        # Find REPCN position in FORMAT
        repcn_pos=$(echo "$format" | tr ':' '\n' | grep -n '^REPCN$' | cut -d':' -f1)
        if [[ -n "$repcn_pos" ]]; then
            repcn=$(echo "$sample" | cut -d':' -f$repcn_pos)
        else
            repcn="N/A"
        fi
        
        # Find GT position
        gt_pos=$(echo "$format" | tr ':' '\n' | grep -n '^GT$' | cut -d':' -f1)
        if [[ -n "$gt_pos" ]]; then
            genotype=$(echo "$sample" | cut -d':' -f$gt_pos)
        else
            genotype="N/A"
        fi
        
        # Extract clinical info from INFO
        disname=$(echo "$line" | grep -oP 'DISNAME=\K[^;]+' || echo "N/A")
        disinher=$(echo "$line" | grep -oP 'DISINHER=\K[^;]+' || echo "N/A")
        disrange=$(echo "$line" | grep -oP 'DISRANGE=\K[^;]+' || echo "N/A")
        
        echo -e "$varid\t$genotype\t$repcn\t$disname\t$disinher\t$disrange"
    done
} > "$SUMMARY_FILE"

echo "  ✓ Summary table created: $SUMMARY_FILE"
echo ""
echo "  Preview:"
column -t -s $'\t' "$SUMMARY_FILE" | head -40

# =============================================================================
# STEP 5: Final Summary
# =============================================================================
echo ""
echo "[Step 5/6] Final Summary"

total=$(grep -v "^#" "$ANNOTATED_VCF" | wc -l)
significant=$(grep -c -E "DISRANGE=.*(Pathogenic|Intermediate)" "$ANNOTATED_VCF" || echo "0")

echo "  ┌─────────────────────────────────────────┐"
echo "  │  SAMPLE: $SAMPLE"
echo "  ├─────────────────────────────────────────┤"
echo "  │  Total loci analyzed:       $total"
echo "  │  Potentially significant:   $significant"
echo "  └─────────────────────────────────────────┘"

if [[ $significant -gt 0 ]]; then
    echo ""
    echo "  ⚠️  CLINICAL ALERT: Potentially significant loci found:"
    grep -E "DISRANGE=.*(Pathogenic|Intermediate)" "$ANNOTATED_VCF" | while read -r line; do
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

if [[ $significant -gt 0 ]]; then
    SIGNIFICANT_ONLY="${EH_PREFIX}_significant_only.tsv"
    
    {
        echo -e "Locus\tGenotype\tRepeatCounts\tDisease\tInheritance\tRange"
        grep -E "DISRANGE=.*(Pathogenic|Intermediate)" "$ANNOTATED_VCF" | while read -r line; do
            varid=$(echo "$line" | grep -oP 'VARID=\K[^;]+' || echo "UNKNOWN")
            
            format=$(echo "$line" | cut -f9)
            sample=$(echo "$line" | cut -f10)
            
            repcn_pos=$(echo "$format" | tr ':' '\n' | grep -n '^REPCN$' | cut -d':' -f1)
            if [[ -n "$repcn_pos" ]]; then
                repcn=$(echo "$sample" | cut -d':' -f$repcn_pos)
            else
                repcn="N/A"
            fi
            
            gt_pos=$(echo "$format" | tr ':' '\n' | grep -n '^GT$' | cut -d':' -f1)
            if [[ -n "$gt_pos" ]]; then
                genotype=$(echo "$sample" | cut -d':' -f$gt_pos)
            else
                genotype="N/A"
            fi
            
            disname=$(echo "$line" | grep -oP 'DISNAME=\K[^;]+' || echo "N/A")
            disinher=$(echo "$line" | grep -oP 'DISINHER=\K[^;]+' || echo "N/A")
            disrange=$(echo "$line" | grep -oP 'DISRANGE=\K[^;]+' || echo "N/A")
            
            echo -e "$varid\t$genotype\t$repcn\t$disname\t$disinher\t$disrange"
        done
    } > "$SIGNIFICANT_ONLY"
    
    echo "  ✓ Significant-only report: $SIGNIFICANT_ONLY"
else
    echo "  ✅ No significant loci - skipping significant-only report"
fi

# =============================================================================
# STEP 6: Plot significant loci
# =============================================================================
echo ""
echo "[Step 6/6] Plotting significant loci"

PLOTS_DIR="${OUTDIR}/reviewer"
mkdir -p "$PLOTS_DIR"

grep -E "DISRANGE=.*(Pathogenic|Intermediate)" "$ANNOTATED_VCF" | grep -oP 'VARID=\K[^;]+' | sort -u | while read -r LOCUS; do
    echo "  → Plotting: $LOCUS"
    
    REViewer \
        --reads "${EH_PREFIX}_realigned.sorted.bam" \
        --vcf "$ANNOTATED_VCF" \
        --reference "$REF" \
        --catalog "$CATALOG" \
        --locus "$LOCUS" \
        --output-prefix "${PLOTS_DIR}/${LOCUS}" \
        2> "${PLOTS_DIR}/reviewer_${LOCUS}.log" || true
    
    # Delete log file if empty
    if [[ ! -s "${PLOTS_DIR}/reviewer_${LOCUS}.log" ]]; then
        rm -f "${PLOTS_DIR}/reviewer_${LOCUS}.log"
    fi
    
    
    # Rename duplicate if needed and move to plots subdirectory
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

echo "  ✓ Plotting complete"

# =============================================================================
# PIPELINE COMPLETE
# =============================================================================
echo ""
echo "================================================================================"
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "================================================================================"
echo "Sample: $SAMPLE"
echo "Summary: $SUMMARY_FILE"
echo "Plots:   $PLOTS_DIR/*.svg"
echo "Completion time: $(date)"
echo "================================================================================"

# List output files
echo ""
echo "Output files:"
ls -lh "$SUMMARY_FILE" 2>/dev/null || echo "  No summary file"
ls -lh "$PLOTS_DIR"/*.svg 2>/dev/null || echo "  No plot files found"
