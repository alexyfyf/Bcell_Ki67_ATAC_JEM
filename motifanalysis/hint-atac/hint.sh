BASE=$(basename $1 .final.bam)
PEAK=$(echo $1 | sed 's/\.final\.bam/_narrow_peaks\.narrowPeak/; s/bams/peaks/')

echo $BASE
echo $PEAK

## call footprint
rgt-hint footprinting --atac-seq --paired-end --organism=mm10 --output-prefix=$BASE $1 $PEAK

## motif matching
rgt-motifanalysis matching --organism=mm10 --input-files ${BASE}.bed

## differential footprinting
# rgt-hint differential --organism=mm10 --bc --nc 12 \
#     --mpbs-files=../match/atac-6mL2_S5_L002_mpbs.bed,../match/atac-6wW2_S2_L001_mpbs.bed \
#     --reads-files=/Myvolume/fyan0011/Lmo2_ATAC/result_atac/shiftedbamfiles/atac-6mL2_S5_L002_shifted_sorted.bam,/Myvolume/fyan0011/Lmo2_ATAC/result_atac/shiftedbamfiles/atac-6wW2_S2_L001_shifted_sorted.bam \
#     --conditions=LSC,WE \
#     --output-location=LSC_WE2

## differential analysis
rgt-hint differential --organism=mm10 --bc --nc 12 \
--mpbs-files=match/A_KO_mpbs.bed,\
match/A_WT_mpbs.bed,\
match/B_KO_mpbs.bed,\
match/B_WT_mpbs.bed,\
match/C_KO_mpbs.bed,\
match/C_WT_mpbs.bed \
--reads-files=../tobias/bams/A_KO.final.bam,../tobias/bams/A_WT.final.bam,\
../tobias/bams/B_KO.final.bam,../tobias/bams/B_WT.final.bam,\
../tobias/bams/C_KO.final.bam,../tobias/bams/C_WT.final.bam \
--conditions=A_KO,A_WT,B_KO,B_WT,C_KO,C_WT \
--output-location=./DiffFootprinting --output-prefix=All