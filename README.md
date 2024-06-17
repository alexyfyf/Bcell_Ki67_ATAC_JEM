# Ki67KO_JEM
ATACseq analysis for the JEM paper

- [Ki67KO\_JEM](#ki67ko_jem)
    - [1, Preprocess using Nextflow pipeline](#1-preprocess-using-nextflow-pipeline)
    - [2, Motif analysis with HINT-ATAC and TOBIAS](#2-motif-analysis-with-hint-atac-and-tobias)
    - [3, Differential analysis with DiffBind](#3-differential-analysis-with-diffbind)
    - [4, Generate tracks with trackhub](#4-generate-tracks-with-trackhub)
    - [Citation](#citation)


### 1, Preprocess using Nextflow pipeline
1. Run custome Nextflow [pipeline](https://github.com/alexyfyf/atac_nf) `preprocess/run_atac_nf.sh` to generate all required files (including BAM, bigwig and narrowPeaks).
2. Software versions contained in [pipeline_info](preprocess/pipeline_info).

### 2, Motif analysis with HINT-ATAC and TOBIAS
1. Generate meta BAM files by merging BAM files (`*.final.bam`) of replicates in each group.
2. Call peaks on meta BAM files for each group using [callpeak.sh](footprint_analysis/tobias/peaks/callpeak.sh).
3. [TOBIAS](https://github.com/loosolab/TOBIAS).
   1. Combine peaks `cat *.narrowPeaks | sort -k1,1 -k2,2n | bedtools merge > merged.bed` as per this [issue](https://github.com/loosolab/TOBIAS/issues/55). 
   2. Run [run_tobias.sh](footprint_analysis/tobias/run_tobias.sh) to prepare corrected signal files and footprint scores.
   3. Run [run_bindetect.sh](footprint_analysis/tobias/run_bindetect.sh) to detect differential transcription factor binding for various comparisons.
4. [HINT-ATAC](https://reg-gen.readthedocs.io/en/latest/hint/introduction.html).
   1. Run [hint.sh](footprint_analysis/hint-atac/hint.sh).

### 3, Differential analysis with DiffBind
1. Run [run_qc.sh](atac_diff/run_qc.sh). 
   1. [samplesheet_atac_ZD.csv](atac_diff/samplesheet_atac_ZD.csv) contains the files been used .
   2. Use [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) to generate consensus peaks and occupancy matrix (`R/3.6.0` and `DiffBind_2.12.0`). 
   3. Plot some QC figures (PCA and clustering).
2. Differential analysis using DESeq2 framework.
   1. Run [1.diff_peak.R](atac_diff/1.diff_peak.R) to get all DMRs and output as RDS, csv and BED files.
3. Run MEME and HOMER on differentially accessible regions (DARs).
   1. Run [run_sizegiven.sh](atac_diff/homer/run_sizegiven.sh) on DMR BED files for HOMER analysis.
   2. Run [convert.sh](atac_diff/meme/convert.sh) on DMR BED files to prepare 500bp fasta files for MEME web analysis.
4. Plotting
   1. [2.plot_v2.R](atac_diff/2.plot_v2.R) to generate figures used in publications (main and supplementary figures).

### 4, Generate tracks with trackhub
1. Put folder `bw` (containing bigwig files) and [generatetrack_bycelltype.py](trackhubs/generatetrack_bycelltype.py) under the same parent directory.
2. Generate trackhub for USCS using [trackhub](https://daler.github.io/trackhub/) python package. Run `python generatetrack_bycelltype.py`
3. Tracks can be viewed at [UCSC track hub](https://genome.ucsc.edu/s/alexyfyf/Ki67_BCell).

### Citation

**Ki67 deficiency impedes chromatin accessibility and BCR gene rearrangement** Zhoujie Ding, Maree Hagan, Feng Yan, Nick Schroer, Jack Polmear, Kim Good-Jacobson, Alexandra Dvorscek, Catherine Pitt, Kristy O'Donnell, Stephen Nutt, Dimitra Zotos, Craig McKenzie, Danika Hill, Marcus Robinson, Isaak Quast, Frank Koentgen, and David Tarlinton. *Journal of Experimental Medicine*. DOI: [10.1084/jem.20232160](https://doi.org/10.1084/jem.20232160).




