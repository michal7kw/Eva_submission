#!/bin/bash

#===============================================================================
# SCRIPT: 11_combine_replicates.sh
# PURPOSE: Combine biological replicates for average signal analysis
#
# DESCRIPTION:
# This script combines biological replicates to create average signal tracks
# and consensus peak sets for improved signal-to-noise ratio and statistical
# power. It merges BAM files, generates average BigWig tracks, creates
# high-confidence consensus peaks, and performs cross-condition signal analysis
# to compare binding patterns across different experimental conditions.
#
# KEY OPERATIONS:
# 1. BAM file merging per condition (TES, TESmut, TEAD1)
# 2. Average BigWig track generation from merged BAMs
# 3. Consensus peak set creation (peaks in ≥2 replicates)
# 4. Signal matrix computation around consensus peaks
# 5. Cross-condition signal analysis and visualization
# 6. Summary statistics generation
#
# METHODOLOGY (Improved for Cut&Tag):
# - Uses samtools for BAM merging and indexing
# - Employs bamCoverage for normalized average signal tracks
# - Applies bedtools multiinter for precise overlap detection
# - Implements IDR (Irreproducible Discovery Rate) for replicate reproducibility
# - Uses summit-based peak consensus for accurate binding sites
# - Utilizes deepTools for signal matrix computation
# - Creates heatmaps and profile plots for visualization
# - Implements stringent consensus criteria with quality weighting
#
# IMPORTANT PARAMETERS:
# - Memory: 16GB (for large BAM file operations)
# - Time: 2 hours (adequate for merging and analysis)
# - Threads: 8 (parallel processing for efficiency)
# - Bin size: 10bp (high-resolution signal tracks)
# - Window: ±3kb around peak centers
# - Consensus threshold: ≥2 replicates
#
# INPUTS:
# - Filtered BAM files from step 4 (per replicate)
# - Individual peak files from step 5 (narrowPeak format)
# - BigWig files from step 6 (per replicate)
# - Three biological replicates per condition
#
# OUTPUTS:
# - Merged BAM files per condition
# - Average BigWig tracks (CPM normalized)
# - Consensus peak sets (high-confidence)
# - Signal matrices around consensus peaks
# - Cross-condition heatmaps and profiles
# - Comprehensive summary statistics
#
# DEPENDENCIES:
# - samtools for BAM manipulation
# - bamCoverage (deepTools) for signal tracks
# - bedtools for genomic interval operations
# - idr package for replicate reproducibility analysis
# - computeMatrix, plotHeatmap (deepTools) for visualization
# - Conda environment: combine_results
#
# USAGE:
# sbatch 11_combine_replicates.sh
#
# NOTES:
# - Improves signal quality through replicate averaging
# - Increases statistical confidence in peak calls
# - Enables robust cross-condition comparisons
# - Essential for publication-quality analysis
# - Provides foundation for downstream comparative studies
#===============================================================================

#SBATCH --job-name=combine_replicates
#SBATCH --output=./logs/11_combine_replicates.out
#SBATCH --error=./logs/11_combine_replicates.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00

echo "=== Combining Replicates for Average Signal Analysis ==="
echo "Start time: $(date)"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate combine_results

# Set up directories
OUTDIR="results/11_combined_replicates_narrow"

# Check required files
REQUIRED_DIRS=("results/04_filtered" "results/05_peaks_narrow" "results/06_bigwig")
for dir in "${REQUIRED_DIRS[@]}"; do
    if [ ! -d "$dir" ]; then
        echo "Error: Required directory $dir not found"
        exit 1
    fi
done

echo "=== Combining BAM files for average signal computation ==="

# Define conditions and replicates (TES and TEAD1 only - TESmut removed)
declare -A CONDITIONS=(
    ["TES"]="TES-1 TES-2 TES-3"
    ["TEAD1"]="TEAD1-1 TEAD1-2 TEAD1-3"
)

# Combine BAM files per condition for average signal
for condition in "${!CONDITIONS[@]}"; do
    echo "Processing $condition condition..."
    
    # Get replicate names
    replicates=(${CONDITIONS[$condition]})
    
    # Check if BAM files exist
    bam_files=""
    for replicate in "${replicates[@]}"; do
        bam_file="results/04_filtered/${replicate}_filtered.bam"
        if [ -f "$bam_file" ]; then
            bam_files="$bam_files $bam_file"
            echo "  Found: $bam_file"
        else
            echo "  Warning: Missing $bam_file"
        fi
    done
    
    if [ -n "$bam_files" ]; then
        echo "  Merging BAM files for $condition..."
        samtools merge -f -@ 8 "$OUTDIR/${condition}_merged.bam" $bam_files
        
        # Index merged BAM
        echo "  Indexing merged BAM for $condition..."
        samtools index "$OUTDIR/${condition}_merged.bam"
        
        # Generate average BigWig from merged BAM
        echo "  Creating average BigWig for $condition..."
        bamCoverage \
            -b "$OUTDIR/${condition}_merged.bam" \
            -o "$OUTDIR/bigwig/${condition}_average.bw" \
            --binSize 10 \
            --normalizeUsing CPM \
            --effectiveGenomeSize 2913022398 \
            --numberOfProcessors 8 \
            --extendReads
            
        echo "  Completed $condition average signal track"
        
        # Calculate read counts for normalization info
        total_reads=$(samtools view -c "$OUTDIR/${condition}_merged.bam")
        echo "  $condition total reads: $total_reads"
        
    else
        echo "  Error: No BAM files found for $condition"
    fi
done

echo "=== Creating consensus peak sets with improved methodology ==="

# Improved consensus peak calling for Cut&Tag data
for condition in "${!CONDITIONS[@]}"; do
    echo "Processing consensus peaks for $condition..."

    # Get individual replicate peak files
    replicates=(${CONDITIONS[$condition]})
    peak_files=()
    replicate_names=()

    for replicate in "${replicates[@]}"; do
        peak_file="results/05_peaks_narrow/${replicate}_peaks.narrowPeak"
        if [ -f "$peak_file" ]; then
            peak_files+=("$peak_file")
            replicate_names+=("$replicate")
            echo "  Found: $peak_file ($(wc -l < $peak_file) peaks)"
        fi
    done

    if [ ${#peak_files[@]} -ge 2 ]; then
        echo "  Creating high-confidence consensus peaks for $condition..."

        # Method 1: Summit-based consensus (most accurate for TF binding)
        # Extract summits from each replicate and expand to ±250bp windows
        summit_files=""
        for i in "${!peak_files[@]}"; do
            # Column 10 in narrowPeak is summit offset from peak start
            awk -v OFS='\t' '{summit=$2+$10; print $1, summit, summit+1, $4, $5, $6, $7, $8, $9}' "${peak_files[$i]}" | \
            bedtools slop -i - -g /beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/hg38.chrom.sizes -b 250 > \
            "$OUTDIR/peaks/${replicate_names[$i]}_summits_expanded.bed"
            summit_files="$summit_files $OUTDIR/peaks/${replicate_names[$i]}_summits_expanded.bed"
        done

        # Use multiinter for precise overlap tracking (only for current condition)
        bedtools multiinter \
            -i $summit_files \
            -names ${replicate_names[@]} \
            > "$OUTDIR/peaks/${condition}_multiinter.bed"

        # High-confidence consensus: present in ≥2 replicates AND ≥50bp overlap
        awk '$4 >= 2 && ($3-$2) >= 50' "$OUTDIR/peaks/${condition}_multiinter.bed" | \
        bedtools merge -i - > "$OUTDIR/peaks/${condition}_consensus_summits.bed"

        # Method 2: Full peak overlap with reciprocal requirement
        # Filter individual peaks first (q-value < 0.05, fold enrichment ≥ 2)
        filtered_files=""
        for i in "${!peak_files[@]}"; do
            awk '$9 >= 1.301 && $7 >= 2' "${peak_files[$i]}" > \
            "$OUTDIR/peaks/${replicate_names[$i]}_filtered.narrowPeak"
            filtered_files="$filtered_files $OUTDIR/peaks/${replicate_names[$i]}_filtered.narrowPeak"
        done

        # Use multiinter on filtered full peaks (only for current condition)
        bedtools multiinter \
            -i $filtered_files \
            -names ${replicate_names[@]} \
            > "$OUTDIR/peaks/${condition}_multiinter_full.bed"

        # Require ≥2 replicates with minimum overlap
        awk '$4 >= 2 && ($3-$2) >= 50' "$OUTDIR/peaks/${condition}_multiinter_full.bed" | \
        bedtools merge -i - -c 4,5 -o max,collapse > "$OUTDIR/peaks/${condition}_consensus_peaks.bed"

        # Method 3: Peak quality-weighted consensus
        # Create scored peaks by averaging signal across overlapping peaks
        for i in "${!peak_files[@]}"; do
            # Extract chr, start, end, score (column 7 = fold enrichment)
            awk -v OFS='\t' '{print $1, $2, $3, $7}' "$OUTDIR/peaks/${replicate_names[$i]}_filtered.narrowPeak" > \
            "$OUTDIR/peaks/${replicate_names[$i]}_scored.bed"
        done

        # Map scores to consensus regions
        bedtools map -a "$OUTDIR/peaks/${condition}_consensus_peaks.bed" \
            -b "$OUTDIR/peaks/${replicate_names[0]}_scored.bed" \
            -c 4 -o mean > "$OUTDIR/peaks/${condition}_consensus_scored.tmp"

        for i in {1..2}; do
            if [ $i -lt ${#replicate_names[@]} ]; then
                bedtools map -a "$OUTDIR/peaks/${condition}_consensus_scored.tmp" \
                    -b "$OUTDIR/peaks/${replicate_names[$i]}_scored.bed" \
                    -c 4 -o mean > "$OUTDIR/peaks/${condition}_consensus_scored.tmp2"
                mv "$OUTDIR/peaks/${condition}_consensus_scored.tmp2" "$OUTDIR/peaks/${condition}_consensus_scored.tmp"
            fi
        done

        # Create final consensus with average scores
        awk -v OFS='\t' '{
            score=0; count=0;
            for(i=4; i<=NF; i++) {
                if($i != ".") {score+=$i; count++}
            }
            avg_score = (count > 0) ? score/count : 0;
            print $1, $2, $3, "peak_"NR, avg_score
        }' "$OUTDIR/peaks/${condition}_consensus_scored.tmp" > "$OUTDIR/peaks/${condition}_consensus_peaks_scored.bed"

        rm "$OUTDIR/peaks/${condition}_consensus_scored.tmp"

        # Copy existing combined peak file for comparison
        if [ -f "results/05_peaks_narrow/${condition}_peaks.narrowPeak" ]; then
            cp "results/05_peaks_narrow/${condition}_peaks.narrowPeak" "$OUTDIR/peaks/${condition}_combined_peaks.narrowPeak"
        fi

        # Report statistics
        summit_consensus=$(wc -l < "$OUTDIR/peaks/${condition}_consensus_summits.bed")
        full_consensus=$(wc -l < "$OUTDIR/peaks/${condition}_consensus_peaks.bed")
        scored_consensus=$(wc -l < "$OUTDIR/peaks/${condition}_consensus_peaks_scored.bed")

        echo "  $condition consensus peaks (summit-based): $summit_consensus"
        echo "  $condition consensus peaks (full overlap): $full_consensus"
        echo "  $condition consensus peaks (quality-scored): $scored_consensus"

    else
        echo "  Warning: Insufficient peak files for $condition consensus (need ≥2, found ${#peak_files[@]})"
    fi
done

echo "=== Running IDR analysis for replicate reproducibility ==="

# IDR (Irreproducible Discovery Rate) analysis for pairwise replicate comparisons
# This is the ENCODE standard for assessing replicate quality
# NOTE: Using relaxed threshold (0.1) and signal.value ranking for Cut&Tag data
mkdir -p "$OUTDIR/idr"

# IDR parameters optimized for Cut&Tag data
IDR_THRESHOLD=0.1  # More relaxed than ChIP-seq default (0.05)
SOFT_IDR_THRESHOLD=0.2  # For rescue ratio calculation

for condition in "${!CONDITIONS[@]}"; do
    echo "Running IDR analysis for $condition..."

    # Get replicate peak files
    replicates=(${CONDITIONS[$condition]})
    peak_files=()

    for replicate in "${replicates[@]}"; do
        peak_file="results/05_peaks_narrow/${replicate}_peaks.narrowPeak"
        if [ -f "$peak_file" ]; then
            peak_files+=("$peak_file")
        fi
    done

    # Run pairwise IDR on all replicate combinations
    if [ ${#peak_files[@]} -eq 3 ]; then
        echo "  Running pairwise IDR comparisons for $condition..."
        echo "  Using IDR threshold: $IDR_THRESHOLD, ranking by signal.value"

        # Rep1 vs Rep2
        idr --samples "${peak_files[0]}" "${peak_files[1]}" \
            --input-file-type narrowPeak \
            --rank signal.value \
            --output-file "$OUTDIR/idr/${condition}_rep1_vs_rep2_idr.txt" \
            --plot \
            --idr-threshold $IDR_THRESHOLD \
            --soft-idr-threshold $SOFT_IDR_THRESHOLD \
            --output-file-type narrowPeak \
            2> "$OUTDIR/idr/${condition}_rep1_vs_rep2_idr.log" || true

        # Rep1 vs Rep3
        idr --samples "${peak_files[0]}" "${peak_files[2]}" \
            --input-file-type narrowPeak \
            --rank signal.value \
            --output-file "$OUTDIR/idr/${condition}_rep1_vs_rep3_idr.txt" \
            --plot \
            --idr-threshold $IDR_THRESHOLD \
            --soft-idr-threshold $SOFT_IDR_THRESHOLD \
            --output-file-type narrowPeak \
            2> "$OUTDIR/idr/${condition}_rep1_vs_rep3_idr.log" || true

        # Rep2 vs Rep3
        idr --samples "${peak_files[1]}" "${peak_files[2]}" \
            --input-file-type narrowPeak \
            --rank signal.value \
            --output-file "$OUTDIR/idr/${condition}_rep2_vs_rep3_idr.txt" \
            --plot \
            --idr-threshold $IDR_THRESHOLD \
            --soft-idr-threshold $SOFT_IDR_THRESHOLD \
            --output-file-type narrowPeak \
            2> "$OUTDIR/idr/${condition}_rep2_vs_rep3_idr.log" || true

        # Check if IDR files exist and have content before processing
        idr_files_exist=true
        for idr_file in "$OUTDIR/idr/${condition}_rep1_vs_rep2_idr.txt" \
                        "$OUTDIR/idr/${condition}_rep1_vs_rep3_idr.txt" \
                        "$OUTDIR/idr/${condition}_rep2_vs_rep3_idr.txt"; do
            if [ ! -f "$idr_file" ] || [ ! -s "$idr_file" ]; then
                idr_files_exist=false
                echo "  Warning: IDR file $idr_file is missing or empty"
            fi
        done

        if [ "$idr_files_exist" = true ]; then
            # Create union of high-confidence IDR peaks
            cat "$OUTDIR/idr/${condition}_rep1_vs_rep2_idr.txt" \
                "$OUTDIR/idr/${condition}_rep1_vs_rep3_idr.txt" \
                "$OUTDIR/idr/${condition}_rep2_vs_rep3_idr.txt" 2>/dev/null | \
            sort -k1,1 -k2,2n | \
            bedtools merge -i - > "$OUTDIR/idr/${condition}_idr_union.bed"

            # Create intersection of high-confidence IDR peaks (most stringent)
            bedtools intersect -a "$OUTDIR/idr/${condition}_rep1_vs_rep2_idr.txt" \
                -b "$OUTDIR/idr/${condition}_rep1_vs_rep3_idr.txt" -u 2>/dev/null | \
            bedtools intersect -a - -b "$OUTDIR/idr/${condition}_rep2_vs_rep3_idr.txt" -u > \
                "$OUTDIR/idr/${condition}_idr_intersection.bed" 2>/dev/null
        else
            echo "  Warning: Skipping IDR union/intersection due to missing files"
            # Create empty placeholder files to prevent downstream errors
            touch "$OUTDIR/idr/${condition}_idr_union.bed"
            touch "$OUTDIR/idr/${condition}_idr_intersection.bed"
        fi

        # Count IDR peaks (handle missing/empty files gracefully)
        idr_12=$(wc -l < "$OUTDIR/idr/${condition}_rep1_vs_rep2_idr.txt" 2>/dev/null || echo 0)
        idr_13=$(wc -l < "$OUTDIR/idr/${condition}_rep1_vs_rep3_idr.txt" 2>/dev/null || echo 0)
        idr_23=$(wc -l < "$OUTDIR/idr/${condition}_rep2_vs_rep3_idr.txt" 2>/dev/null || echo 0)
        idr_union=$(wc -l < "$OUTDIR/idr/${condition}_idr_union.bed" 2>/dev/null || echo 0)
        idr_intersect=$(wc -l < "$OUTDIR/idr/${condition}_idr_intersection.bed" 2>/dev/null || echo 0)

        echo "  $condition IDR results:"
        echo "    Rep1 vs Rep2 (IDR<$IDR_THRESHOLD): $idr_12 peaks"
        echo "    Rep1 vs Rep3 (IDR<$IDR_THRESHOLD): $idr_13 peaks"
        echo "    Rep2 vs Rep3 (IDR<$IDR_THRESHOLD): $idr_23 peaks"
        echo "    Union of all IDR peaks: $idr_union peaks"
        echo "    Intersection (most stringent): $idr_intersect peaks"

        # Check IDR logs for diagnostic info
        if [ $idr_12 -eq 0 ] && [ $idr_13 -eq 0 ] && [ $idr_23 -eq 0 ]; then
            echo "  Note: IDR returned 0 peaks - checking logs for details..."
            grep -i "warning\|error\|peaks\|number" "$OUTDIR/idr/${condition}_rep1_vs_rep2_idr.log" 2>/dev/null | head -5 || true
        fi

    elif [ ${#peak_files[@]} -eq 2 ]; then
        echo "  Running IDR for 2 replicates of $condition..."
        idr --samples "${peak_files[0]}" "${peak_files[1]}" \
            --input-file-type narrowPeak \
            --rank signal.value \
            --output-file "$OUTDIR/idr/${condition}_idr.txt" \
            --plot \
            --idr-threshold $IDR_THRESHOLD \
            --soft-idr-threshold $SOFT_IDR_THRESHOLD \
            --output-file-type narrowPeak \
            2> "$OUTDIR/idr/${condition}_idr.log" || true

        idr_count=$(wc -l < "$OUTDIR/idr/${condition}_idr.txt" 2>/dev/null || echo 0)
        echo "  $condition IDR peaks (IDR<$IDR_THRESHOLD): $idr_count"
    else
        echo "  Warning: Need at least 2 replicates for IDR analysis"
    fi
done

echo "=== Creating signal matrices around consensus peaks ==="

# Create signal matrices for each condition using their own consensus peaks
# Use the quality-scored consensus peaks for best results
for condition in "${!CONDITIONS[@]}"; do
    if [ -f "$OUTDIR/peaks/${condition}_consensus_peaks_scored.bed" ] && [ -f "$OUTDIR/bigwig/${condition}_average.bw" ]; then
        echo "Computing signal matrix for $condition around its consensus peaks..."

        computeMatrix reference-point \
            --referencePoint center \
            -b 3000 -a 3000 \
            -R "$OUTDIR/peaks/${condition}_consensus_peaks_scored.bed" \
            -S "$OUTDIR/bigwig/${condition}_average.bw" \
            --skipZeros \
            -o "$OUTDIR/matrices/${condition}_self_matrix.gz" \
            --outFileSortedRegions "$OUTDIR/matrices/${condition}_self_regions.bed" \
            -p 8
    fi
done

echo "=== Cross-condition signal analysis ==="

# Create matrices showing each condition's signal on all other conditions' peaks
conditions=("TES" "TEAD1")  # TESmut removed - failed sample
for peak_condition in "${conditions[@]}"; do
    if [ -f "$OUTDIR/peaks/${peak_condition}_consensus_peaks.bed" ]; then
        
        # Collect all BigWig files
        bigwig_files=""
        bigwig_labels=""
        
        for signal_condition in "${conditions[@]}"; do
            if [ -f "$OUTDIR/bigwig/${signal_condition}_average.bw" ]; then
                bigwig_files="$bigwig_files $OUTDIR/bigwig/${signal_condition}_average.bw"
                bigwig_labels="$bigwig_labels ${signal_condition}_avg"
            fi
        done
        
        if [ -n "$bigwig_files" ]; then
            echo "Computing cross-condition matrix for ${peak_condition} peaks..."
            
            computeMatrix reference-point \
                --referencePoint center \
                -b 3000 -a 3000 \
                -R "$OUTDIR/peaks/${peak_condition}_consensus_peaks.bed" \
                -S $bigwig_files \
                --samplesLabel $bigwig_labels \
                --skipZeros \
                -o "$OUTDIR/matrices/${peak_condition}_peaks_all_signals.gz" \
                --outFileSortedRegions "$OUTDIR/matrices/${peak_condition}_peaks_all_regions.bed" \
                -p 8
                
            # Create heatmap visualization
            plotHeatmap \
                -m "$OUTDIR/matrices/${peak_condition}_peaks_all_signals.gz" \
                -out "$OUTDIR/${peak_condition}_peaks_signal_heatmap.pdf" \
                --plotTitle "${peak_condition} peaks: Signal from all conditions" \
                --whatToShow 'heatmap and colorbar' \
                --colorMap RdBu_r \
                --zMin -2 --zMax 2
                
            # Create profile plot
            plotProfile \
                -m "$OUTDIR/matrices/${peak_condition}_peaks_all_signals.gz" \
                -out "$OUTDIR/${peak_condition}_peaks_signal_profile.pdf" \
                --plotTitle "${peak_condition} peaks: Average signal profiles" \
                --plotType lines \
                --perGroup
        fi
    fi
done

echo "=== Generating comprehensive summary statistics ==="

# Create summary file
SUMMARY_FILE="$OUTDIR/REPLICATE_COMBINATION_SUMMARY_NARROW.txt"
echo "Cut&Tag Replicate Combination Analysis Summary (NARROW PEAKS - IMPROVED)" > "$SUMMARY_FILE"
echo "Generated on: $(date)" >> "$SUMMARY_FILE"
echo "=========================================================================" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

for condition in "${!CONDITIONS[@]}"; do
    echo "=== $condition Condition ===" >> "$SUMMARY_FILE"

    # Read counts from merged BAM
    if [ -f "$OUTDIR/${condition}_merged.bam" ]; then
        reads=$(samtools view -c "$OUTDIR/${condition}_merged.bam")
        echo "Total reads in merged BAM: $reads" >> "$SUMMARY_FILE"
    fi

    echo "" >> "$SUMMARY_FILE"
    echo "Peak Calling Results:" >> "$SUMMARY_FILE"
    echo "---------------------" >> "$SUMMARY_FILE"

    # Individual replicate peak counts
    replicates=(${CONDITIONS[$condition]})
    for replicate in "${replicates[@]}"; do
        peak_file="results/05_peaks_narrow/${replicate}_peaks.narrowPeak"
        if [ -f "$peak_file" ]; then
            count=$(wc -l < "$peak_file")
            echo "  $replicate: $count peaks" >> "$SUMMARY_FILE"
        fi
    done

    echo "" >> "$SUMMARY_FILE"
    echo "Consensus Peak Sets:" >> "$SUMMARY_FILE"
    echo "--------------------" >> "$SUMMARY_FILE"

    # Summit-based consensus
    if [ -f "$OUTDIR/peaks/${condition}_consensus_summits.bed" ]; then
        summit_peaks=$(wc -l < "$OUTDIR/peaks/${condition}_consensus_summits.bed")
        echo "  Summit-based consensus (±250bp): $summit_peaks peaks" >> "$SUMMARY_FILE"
    fi

    # Full peak consensus
    if [ -f "$OUTDIR/peaks/${condition}_consensus_peaks.bed" ]; then
        consensus_peaks=$(wc -l < "$OUTDIR/peaks/${condition}_consensus_peaks.bed")
        echo "  Full overlap consensus (≥2 reps, ≥50bp): $consensus_peaks peaks" >> "$SUMMARY_FILE"
    fi

    # Quality-scored consensus
    if [ -f "$OUTDIR/peaks/${condition}_consensus_peaks_scored.bed" ]; then
        scored_peaks=$(wc -l < "$OUTDIR/peaks/${condition}_consensus_peaks_scored.bed")
        avg_score=$(awk '{sum+=$5; count++} END {print sum/count}' "$OUTDIR/peaks/${condition}_consensus_peaks_scored.bed")
        echo "  Quality-scored consensus: $scored_peaks peaks (avg score: $avg_score)" >> "$SUMMARY_FILE"
    fi

    # Combined peaks (original method for comparison)
    if [ -f "$OUTDIR/peaks/${condition}_combined_peaks.narrowPeak" ]; then
        combined_peaks=$(wc -l < "$OUTDIR/peaks/${condition}_combined_peaks.narrowPeak")
        echo "  Combined peaks (original MACS2): $combined_peaks peaks" >> "$SUMMARY_FILE"
    fi

    echo "" >> "$SUMMARY_FILE"
    echo "IDR Reproducibility Analysis:" >> "$SUMMARY_FILE"
    echo "------------------------------" >> "$SUMMARY_FILE"

    # IDR results
    if [ -f "$OUTDIR/idr/${condition}_rep1_vs_rep2_idr.txt" ]; then
        idr_12=$(wc -l < "$OUTDIR/idr/${condition}_rep1_vs_rep2_idr.txt" 2>/dev/null || echo 0)
        idr_13=$(wc -l < "$OUTDIR/idr/${condition}_rep1_vs_rep3_idr.txt" 2>/dev/null || echo 0)
        idr_23=$(wc -l < "$OUTDIR/idr/${condition}_rep2_vs_rep3_idr.txt" 2>/dev/null || echo 0)

        echo "  Rep1 vs Rep2 (IDR<0.1): $idr_12 peaks" >> "$SUMMARY_FILE"
        echo "  Rep1 vs Rep3 (IDR<0.1): $idr_13 peaks" >> "$SUMMARY_FILE"
        echo "  Rep2 vs Rep3 (IDR<0.1): $idr_23 peaks" >> "$SUMMARY_FILE"

        # Average reproducibility
        avg_idr=$(echo "scale=0; ($idr_12 + $idr_13 + $idr_23) / 3" | bc)
        echo "  Average pairwise IDR peaks: $avg_idr" >> "$SUMMARY_FILE"

        if [ -f "$OUTDIR/idr/${condition}_idr_union.bed" ]; then
            union=$(wc -l < "$OUTDIR/idr/${condition}_idr_union.bed")
            echo "  Union of all IDR peaks: $union peaks" >> "$SUMMARY_FILE"
        fi

        if [ -f "$OUTDIR/idr/${condition}_idr_intersection.bed" ]; then
            intersect=$(wc -l < "$OUTDIR/idr/${condition}_idr_intersection.bed")
            echo "  Intersection (most stringent): $intersect peaks" >> "$SUMMARY_FILE"
        fi
    elif [ -f "$OUTDIR/idr/${condition}_idr.txt" ]; then
        idr_count=$(wc -l < "$OUTDIR/idr/${condition}_idr.txt")
        echo "  IDR peaks (2 replicates, IDR<0.1): $idr_count" >> "$SUMMARY_FILE"
    else
        echo "  IDR analysis not performed (insufficient replicates)" >> "$SUMMARY_FILE"
    fi

    echo "" >> "$SUMMARY_FILE"
    echo "Quality Metrics:" >> "$SUMMARY_FILE"
    echo "----------------" >> "$SUMMARY_FILE"

    # Calculate reproducibility rate (IDR peaks / average individual peaks)
    if [ -f "$OUTDIR/idr/${condition}_rep1_vs_rep2_idr.txt" ]; then
        rep1_peaks=$(wc -l < "results/05_peaks_narrow/${replicates[0]}_peaks.narrowPeak" 2>/dev/null || echo 0)
        rep2_peaks=$(wc -l < "results/05_peaks_narrow/${replicates[1]}_peaks.narrowPeak" 2>/dev/null || echo 0)
        rep3_peaks=$(wc -l < "results/05_peaks_narrow/${replicates[2]}_peaks.narrowPeak" 2>/dev/null || echo 0)

        if [ $rep1_peaks -gt 0 ] && [ $rep2_peaks -gt 0 ]; then
            avg_peaks=$(echo "scale=0; ($rep1_peaks + $rep2_peaks + $rep3_peaks) / 3" | bc)
            repro_rate=$(echo "scale=2; ($avg_idr / $avg_peaks) * 100" | bc)
            echo "  Reproducibility rate: ${repro_rate}% (avg IDR / avg individual)" >> "$SUMMARY_FILE"
        fi
    fi

    echo "" >> "$SUMMARY_FILE"
done

echo "" >> "$SUMMARY_FILE"
echo "Recommended Peak Sets for Downstream Analysis:" >> "$SUMMARY_FILE"
echo "-----------------------------------------------" >> "$SUMMARY_FILE"
echo "1. MOST STRINGENT: IDR intersection peaks (highest confidence)" >> "$SUMMARY_FILE"
echo "   Files: idr/*_idr_intersection.bed" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
echo "2. RECOMMENDED: Quality-scored consensus peaks (balanced stringency)" >> "$SUMMARY_FILE"
echo "   Files: peaks/*_consensus_peaks_scored.bed" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
echo "3. PERMISSIVE: IDR union or summit-based consensus (maximum coverage)" >> "$SUMMARY_FILE"
echo "   Files: idr/*_idr_union.bed or peaks/*_consensus_summits.bed" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

echo "Files created:"
echo "- Average BigWig tracks: $OUTDIR/bigwig/"
echo "- Consensus peak sets: $OUTDIR/peaks/"
echo "- Signal matrices: $OUTDIR/matrices/"
echo "- Cross-condition heatmaps: $OUTDIR/*_heatmap.pdf"
echo "- Summary report: $SUMMARY_FILE"

echo "=== Narrow peak replicate combination completed ==="
echo "End time: $(date)"