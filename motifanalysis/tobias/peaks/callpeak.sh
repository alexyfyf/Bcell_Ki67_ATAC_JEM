bam=$1
echo $bam

name=$(basename $bam .final.bam)
echo $name 

species='mm'
bedtools bamtobed -i $bam > ${name}_pe.bed 
macs2 callpeak -t ${name}_pe.bed -n ${name}_narrow -f BED -g ${species} -q 0.01 --nomodel --shift -75 --extsize 150 --call-summits --keep-dup all

