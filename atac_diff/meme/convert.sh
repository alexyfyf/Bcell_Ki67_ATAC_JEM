awk -vOFS="\t" 'width=$3-$2 {if(width % 2 != 0) {width+=1} ; mid=$2+width/2; print $1,mid,mid+1,$4}' $1 > $(basename $1 .bed)_summits.bed

bedtools slop -b 250 -i $(basename $1 .bed)_summits.bed -g mm10.chrom.sizes > $(basename $1 .bed)_500bp.bed

bedtools getfasta -fi /home/fyan0011/ls25/feng.yan_raw/refFiles/ucsc_mm10/mm10.fa -bed $(basename $1 .bed)_500bp.bed > $(basename $1 .bed)_500bp.fasta
