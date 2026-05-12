#!/usr/bin/env bash
# =============================================================================
# Short-read Multi-locus STR Analysis Pipeline
# ExpansionHunter + STRipy + REViewer
# =============================================================================
# Purpose: Detect pathogenic repeat expansions from Illumina short-read WGS data
#          (Multi-locus analysis - analyzes all loci in the catalog)
# =============================================================================
# Input:
#   - BAM/CRAM file (aligned to hg38)
#   - Reference FASTA (hg38)
#   - Variant catalog JSON (ExpansionHunter format)
# =============================================================================
# Output:
#   - ${SAMPLE}_multi.vcf                    → Raw ExpansionHunter output
#   - ${SAMPLE}_multi_annotated.vcf         → STRipy annotated VCF
#   - ${SAMPLE}_multi_info.tsv              → Full table (ALL loci)
#   - ${SAMPLE}_significant_only_info.tsv   → Only Pathogenic/Intermediate loci
#   - ${SAMPLE}_multi_realigned.sorted.bam  → Realigned BAMlet for REViewer
#   - reviewer/${SAMPLE}_*.svg              → REViewer plots for significant loci
# =============================================================================
# Dependencies: ExpansionHunter, REViewer, samtools, curl
# =============================================================================

set -euo pipefail

# =============================================================================
# Default values
# =============================================================================
SAMPLE=""
INPUT_BAM=""
REF=""
CATALOG=""
OUTDIR=""
THREADS=12

# =============================================================================
# Usage
# =============================================================================
usage() {
    cat << EOF
Usage: $(basename "$0") --input <BAM> --ref <FASTA> --catalog <JSON> --output <DIR> [OPTIONS]

Required:
  --input     Input BAM/CRAM file
  --ref       Reference genome FASTA
  --catalog   Variant catalog JSON
  --output    Output directory

Optional:
  --sample    Sample ID
  --threads   Number of threads (default: 12)
EOF
    exit 1
}

# =============================================================================
# Argument parsing
# =============================================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)   INPUT_BAM="$2"; shift 2 ;;
        --ref)     REF="$2";       shift 2 ;;
        --catalog) CATALOG="$2";   shift 2 ;;
        --output)  OUTDIR="$2";    shift 2 ;;
        --sample)  SAMPLE="$2";    shift 2 ;;
        --threads) THREADS="$2";   shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# Validation
if [[ -z "$INPUT_BAM" || -z "$REF" || -z "$CATALOG" || -z "$OUTDIR" ]]; then
    echo "❌ Error: Missing required arguments."
    usage
fi

[[ ! -f "$INPUT_BAM" ]] && { echo "❌ BAM not found"; exit 1; }
[[ ! -f "$REF" ]]       && { echo "❌ Reference not found"; exit 1; }
[[ ! -f "$CATALOG" ]]   && { echo "❌ Catalog not found"; exit 1; }

if [[ -z "$SAMPLE" ]]; then
    SAMPLE=$(basename "$INPUT_BAM" | sed 's/\.[^.]*$//')
fi

# =============================================================================
# Setup
# =============================================================================
mkdir -p "$OUTDIR"/reviewer
cd "$OUTDIR"

EH_PREFIX="${SAMPLE}_multi"
ANNOTATED_VCF="${EH_PREFIX}_annotated.vcf"
SUMMARY_FILE="${SAMPLE}_multi_info.tsv"
SIG_REPORT="${SAMPLE}_significant_only_info.tsv"
PLOTS_DIR="reviewer"

echo "================================================================================"
echo "STR Expansion Analysis Pipeline"
echo "Sample   : $SAMPLE"
echo "Started  : $(date)"
echo "================================================================================"

# =============================================================================
# 1. ExpansionHunter
# =============================================================================
echo ""
echo "[1/6] Running ExpansionHunter..."
if [[ -s "${EH_PREFIX}.vcf" ]]; then
    echo " ✓ VCF already exists → skipping"
else
    ExpansionHunter --reads "$INPUT_BAM" --reference "$REF" --variant-catalog "$CATALOG" \
        --output-prefix "$EH_PREFIX" --threads "$THREADS" --analysis-mode streaming
    echo " ✓ ExpansionHunter completed"
fi

# =============================================================================
# 2. Prepare BAMlet
# =============================================================================
echo ""
echo "[2/6] Preparing BAMlet for REViewer..."
if [[ ! -f "${EH_PREFIX}_realigned.sorted.bam" ]]; then
    samtools sort -@ "$((THREADS/2))" -o "${EH_PREFIX}_realigned.sorted.bam" "${EH_PREFIX}_realigned.bam"
    samtools index "${EH_PREFIX}_realigned.sorted.bam"
    echo " ✓ BAMlet ready"
else
    echo " ✓ BAMlet already exists"
fi

# =============================================================================
# 3. STRipy Annotation
# =============================================================================
echo ""
echo "[3/6] Annotating with STRipy..."
if [[ ! -s "$ANNOTATED_VCF" ]]; then
    curl -s -F "file=@${EH_PREFIX}.vcf" https://api.stripy.org/annotateVCF > "$ANNOTATED_VCF"
    echo " ✓ Annotation completed"
else
    echo " ✓ Annotated VCF already exists"
fi

# =============================================================================
# 4. Full Summary Table - CORRECT TAB-SEPARATED COLUMNS
# =============================================================================
echo ""
echo "[4/6] Creating full summary table..."

# Create header with proper TAB separation
printf "#CHROM\tPOS\tFILTER\tEND\tReference_Copy_Number\tReference_Length_bp\tRepeat_Unit\tVariant_ID\tDisease_Name\tInheritance\tClinical_Range\tGenotype\tSupport_Type\tRepeat_Count\tConfidence_Interval\tSpanning_Reads\tFlanking_Reads\tInRepeat_Reads\tLocus_Coverage\n" > "$SUMMARY_FILE"

awk -F'\t' '!/^#/ && NF>=10 {
    if($1=="" || $2=="") next;

    varid=disname=disinher=disrange=end=ref_cn=rl=ru="";
    split($8, info, ";");
    for(i in info) {
        eq=index(info[i],"=");
        if(eq>0) {
            k=substr(info[i],1,eq-1);
            v=substr(info[i],eq+1);
            if(k=="END") end=v;
            else if(k=="REF") ref_cn=v;
            else if(k=="RL") rl=v;
            else if(k=="RU") ru=v;
            else if(k=="VARID") varid=v;
            else if(k=="DISNAME") disname=v;
            else if(k=="DISINHER") disinher=v;
            else if(k=="DISRANGE") disrange=v;
        }
    }
    if(varid=="") next;

    split($9,fmt,":"); split($10,samp,":");
    gt=so=repcn=repci=adsp=adfl=adir=lc="./.";
    for(i in fmt) {
        val=samp[i];
        if(fmt[i]=="GT") gt=val;
        else if(fmt[i]=="SO") so=val;
        else if(fmt[i]=="REPCN") repcn=val;
        else if(fmt[i]=="REPCI") repci=val;
        else if(fmt[i]=="ADSP") adsp=val;
        else if(fmt[i]=="ADFL") adfl=val;
        else if(fmt[i]=="ADIR") adir=val;
        else if(fmt[i]=="LC") lc=val;
    }

    print $1 "\t" $2 "\t" $7 "\t" end "\t" ref_cn "\t" rl "\t" ru "\t" varid "\t" disname "\t" disinher "\t" disrange "\t" gt "\t" so "\t" repcn "\t" repci "\t" adsp "\t" adfl "\t" adir "\t" lc;
}' "$ANNOTATED_VCF" >> "$SUMMARY_FILE"

echo " ✓ Full summary table created: $SUMMARY_FILE"

# =============================================================================
# 5. Significant Only Report
# =============================================================================
echo ""
echo "[5/6] Creating significant-only report..."

significant=$(grep -c -E "DISRANGE=.*(Pathogenic|Intermediate)" "$ANNOTATED_VCF" 2>/dev/null || echo 0)

if [[ $significant -gt 0 ]]; then
    # Correct header with proper TAB separation
    printf "#CHROM\tPOS\tFILTER\tEND\tReference_Copy_Number\tReference_Length_bp\tRepeat_Unit\tVariant_ID\tDisease_Name\tInheritance\tClinical_Range\tGenotype\tSupport_Type\tRepeat_Count\tConfidence_Interval\tSpanning_Reads\tFlanking_Reads\tInRepeat_Reads\tLocus_Coverage\n" > "$SIG_REPORT"

    awk -F'\t' '!/^#/ && $8 ~ /Pathogenic|Intermediate/ {
        varid=disname=disinher=disrange=end=ref_cn=rl=ru="";
        split($8, info, ";");
        for(i in info) {
            eq=index(info[i],"=");
            if(eq>0) {
                k=substr(info[i],1,eq-1);
                v=substr(info[i],eq+1);
                if(k=="END") end=v;
                else if(k=="REF") ref_cn=v;
                else if(k=="RL") rl=v;
                else if(k=="RU") ru=v;
                else if(k=="VARID") varid=v;
                else if(k=="DISNAME") disname=v;
                else if(k=="DISINHER") disinher=v;
                else if(k=="DISRANGE") disrange=v;
            }
        }

        split($9,fmt,":"); split($10,samp,":");
        gt=so=repcn=repci=adsp=adfl=adir=lc="./.";
        for(i in fmt) {
            val=samp[i];
            if(fmt[i]=="GT") gt=val;
            else if(fmt[i]=="SO") so=val;
            else if(fmt[i]=="REPCN") repcn=val;
            else if(fmt[i]=="REPCI") repci=val;
            else if(fmt[i]=="ADSP") adsp=val;
            else if(fmt[i]=="ADFL") adfl=val;
            else if(fmt[i]=="ADIR") adir=val;
            else if(fmt[i]=="LC") lc=val;
        }

        print $1 "\t" $2 "\t" $7 "\t" end "\t" ref_cn "\t" rl "\t" ru "\t" varid "\t" disname "\t" disinher "\t" disrange "\t" gt "\t" so "\t" repcn "\t" repci "\t" adsp "\t" adfl "\t" adir "\t" lc;
    }' "$ANNOTATED_VCF" >> "$SIG_REPORT"

    echo " ✓ Significant-only report created: $SIG_REPORT ($significant loci)"
else
    echo " ✅ No significant loci found - skipping significant-only report"
fi

# =============================================================================
# 6. REViewer Plots
# =============================================================================
echo ""
echo "[6/6] Generating REViewer plots..."
mkdir -p "$PLOTS_DIR"

if [[ $significant -gt 0 ]]; then
    grep -E "DISRANGE=.*(Pathogenic|Intermediate)" "$ANNOTATED_VCF" | grep -oP 'VARID=\K[^;]+' | sort -u | while read -r LOCUS; do
        echo " → Plotting $LOCUS"
        REViewer \
            --reads "${EH_PREFIX}_realigned.sorted.bam" \
            --vcf "$ANNOTATED_VCF" \
            --reference "$REF" \
            --catalog "$CATALOG" \
            --locus "$LOCUS" \
            --output-prefix "${PLOTS_DIR}/${LOCUS}" || true
    done
fi

# =============================================================================
# Final Summary
# =============================================================================
echo ""
echo "================================================================================"
echo "✅ PIPELINE COMPLETED SUCCESSFULLY"
echo "Sample               : $SAMPLE"
echo "Significant findings : $significant"
echo "Output folder        : $OUTDIR"
echo "================================================================================"

exit 0
