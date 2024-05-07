for i in {A..C}
do 
TOBIAS BINDetect --motifs /Myvolume/fyan0011/Lmo2_ATAC/tobias_result/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt --signals ${i}_KO/${i}_KO_footprints.bw ${i}_WT/${i}_WT_footprints.bw --genome /Myvolume/fyan0011/ref/mmu/ucsc_mm10/mm10.fa --peaks peaks/merged_clipped.bed --outdir ${i}_KO2WT_output --cores 12 --cond-names ${i}_KO ${i}_WT &> ${i}_bindetect.log
done

## time series
TOBIAS BINDetect --motifs /Myvolume/fyan0011/Lmo2_ATAC/tobias_result/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt --signals A_WT/A_WT_footprints.bw B_WT/B_WT_footprints.bw C_WT/C_WT_footprints.bw --genome /Myvolume/fyan0011/ref/mmu/ucsc_mm10/mm10.fa --peaks peaks/merged_clipped.bed --outdir WT_ABC_output --time-series --cores 12 --cond-names A B C &> WT_bindetect.log

TOBIAS BINDetect --motifs /Myvolume/fyan0011/Lmo2_ATAC/tobias_result/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt --signals A_KO/A_KO_footprints.bw B_KO/B_KO_footprints.bw C_KO/C_KO_footprints.bw --genome /Myvolume/fyan0011/ref/mmu/ucsc_mm10/mm10.fa --peaks peaks/merged_clipped.bed --outdir KO_ABC_output --time-series --cores 12 --cond-names A B C &> KO_bindetect.log 
