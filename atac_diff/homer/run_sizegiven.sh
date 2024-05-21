awk -v OFS='\t' '{print NR, $1, $2, $3, $6}'  $1 > $(basename $1)

findMotifsGenome.pl $(basename $1) /home/fyan0011/ls25/feng.yan_raw/refFiles/ucsc_mm10/mm10.fa $(basename $1 .bed) -size given -len 8 -p 12
