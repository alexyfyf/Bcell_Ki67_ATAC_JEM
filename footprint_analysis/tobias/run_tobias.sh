BASE=$(basename $1 .final.bam)
echo $BASE

TOBIAS ATACorrect --bam $1 --genome /Myvolume/fyan0011/ref/mmu/ucsc_mm10/mm10.fa --peaks peaks/merged.bed --outdir $BASE --cores 12

TOBIAS FootprintScores --signal ${BASE}/${BASE}.final_corrected.bw --regions  peaks/merged.bed --output ${BASE}/${BASE}_footprints.bw --cores 12
