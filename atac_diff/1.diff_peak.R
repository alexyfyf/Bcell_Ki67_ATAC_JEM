
db <- readRDS("data/db.rds")
db$samples$Treatment <- db$samples$Tissue

## use function to plot PCA and hclust
source("pca_clustering.R")

pca_clustering(db, sampleid = 1:24, prefix = "top5k", top = 5000)
pca_clustering(db, sampleid = 1:24, prefix = "top10k", top = 10000)
pca_clustering(db, sampleid = 1:24, prefix = "top1k", top = 1000)

pca_clustering(db, sampleid = 1:8, prefix = "ppb_5k", top = 5000)
pca_clustering(db, sampleid = 9:16, prefix = "pb_5k", top = 5000)
pca_clustering(db, sampleid = 17:24, prefix = "spb_5k", top = 5000)

peaks <- db$peaks[[1]] %>% makeGRangesFromDataFrame()
export.bed(peaks, 'output/consensus_peaks.bed')

## consensus peak distribution
source("comp_meth_anno.R")
comp_meth_anno(peaks %>% mutate(Fold=1), cpgs_info = NULL, simplify = T, column = "Fold")
ggsave("output/consensus_anno_pie.pdf", width = 5, height = 7)

comp_meth_anno(peaks %>% mutate(Fold=1), cpgs_info = NULL, simplify = T, column = "Fold",plot = "bar")
ggsave("output/consensus_anno_bar.pdf", width = 7, height = 7)

## KO vs WT
## A: pre-pro-B
db.ppb = dba(db, mask = db$masks$`Pre-pro-B`)
dar.ppb = dba.contrast(db.ppb, group1 = db.ppb$masks$Ki67KO,
                     group2 = db.ppb$masks$WT) %>%
  dba.analyze(., method = DBA_ALL_METHODS)

dba.show(dar.ppb, bContrasts = T)

dar.ppb.report=dba.report(dar.ppb)

write.csv(data.frame(dar.ppb.report),"output/dar.ppb.csv")

export.bed(dar.ppb.report, 'output/dar.ppb.bed')
export.bed(dar.ppb.report %>% filter(Fold>0), 'output/dar.ppb_up.bed')
export.bed(dar.ppb.report %>% filter(Fold<0), 'output/dar.ppb_down.bed')

saveRDS(dar.ppb.report,"output/dar.ppb.report.rds")

dar.ppb.anno <- comp_meth_anno(dar.ppb.report, cpgs_info = NULL, column = "Fold")
ggsave("output/dar.ppb_annotation_pie.pdf", width = 10, height = 7)

## B: pro-B
db.pb = dba(db, mask = db$masks$`Pro-B`)
dar.pb = dba.contrast(db.pb, group1 = db.pb$masks$Ki67KO,
                       group2 = db.pb$masks$WT) %>%
  dba.analyze(., method = DBA_ALL_METHODS)

dba.show(dar.pb, bContrasts = T)

dar.pb.report=dba.report(dar.pb)

write.csv(data.frame(dar.pb.report),"output/dar.pb.csv")

export.bed(dar.pb.report, 'output/dar.pb.bed')
export.bed(dar.pb.report %>% filter(Fold>0), 'output/dar.pb_up.bed')
export.bed(dar.pb.report %>% filter(Fold<0), 'output/dar.pb_down.bed')

saveRDS(dar.pb.report,"output/dar.pb.report.rds")

dar.pb.anno <- comp_meth_anno(dar.pb.report, cpgs_info = NULL, column = "Fold")
ggsave("output/dar.pb_annotation_pie.pdf", width = 10, height = 7)

## C: small pre-B
db.spb = dba(db, mask = db$masks$`small pre-B`)
dar.spb = dba.contrast(db.spb, group1 = db.spb$masks$Ki67KO,
                       group2 = db.spb$masks$WT) %>%
  dba.analyze(., method = DBA_ALL_METHODS)

dba.show(dar.spb, bContrasts = T)

dar.spb.report=dba.report(dar.spb)

write.csv(data.frame(dar.spb.report),"output/dar.spb.csv")

export.bed(dar.spb.report, 'output/dar.spb.bed')
export.bed(dar.spb.report %>% filter(Fold>0), 'output/dar.spb_up.bed')
export.bed(dar.spb.report %>% filter(Fold<0), 'output/dar.spb_down.bed')

saveRDS(dar.spb.report,"output/dar.spb.report.rds")

dar.spb.anno <- comp_meth_anno(dar.spb.report, cpgs_info = NULL, column = "Fold")
ggsave("output/dar.spb_annotation_pie.pdf", width = 10, height = 7)

## biological comparison

## pro vs pre pro, small pre vs pro

## wt1 proB vs pre-pro B
db.wt = dba(db, mask = db$masks$WT)
db.wt1 = dba.contrast(db.wt, group1 = db.wt$masks$`Pro-B`,
                       group2 = db.wt$masks$`Pre-pro-B`) %>%
  dba.analyze(., method = DBA_ALL_METHODS)
## wt2 small pre B vs proB
db.wt2 = dba.contrast(db.wt, group1 = db.wt$masks$`small pre-B`,
                      group2 = db.wt$masks$`Pro-B`) %>%
  dba.analyze(., method = DBA_ALL_METHODS)

dba.show(db.wt1, bContrasts = T)
dba.show(db.wt2, bContrasts = T)

dar.wt1.report=dba.report(db.wt1)

write.csv(data.frame(dar.wt1.report),"output/dar.wt1.csv")

export.bed(dar.wt1.report, 'output/dar.wt1.bed')
export.bed(dar.wt1.report %>% filter(Fold>0), 'output/dar.wt1_up.bed')
export.bed(dar.wt1.report %>% filter(Fold<0), 'output/dar.wt1_down.bed')

saveRDS(dar.wt1.report,"output/dar.wt1.report.rds")

dar.wt1.anno <- comp_meth_anno(dar.wt1.report, cpgs_info = NULL, column = "Fold")
ggsave("output/dar.wt1_annotation_pie.pdf", width = 10, height = 7)

dar.wt2.report=dba.report(db.wt2)

write.csv(data.frame(dar.wt2.report),"output/dar.wt2.csv")

export.bed(dar.wt2.report, 'output/dar.wt2.bed')
export.bed(dar.wt2.report %>% filter(Fold>0), 'output/dar.wt2_up.bed')
export.bed(dar.wt2.report %>% filter(Fold<0), 'output/dar.wt2_down.bed')

saveRDS(dar.wt2.report,"output/dar.wt2.report.rds")

dar.wt2.anno <- comp_meth_anno(dar.wt2.report, cpgs_info = NULL, column = "Fold")
ggsave("output/dar.wt2_annotation_pie.pdf", width = 10, height = 7)

db.ko = dba(db, mask = db$masks$Ki67KO)
## ko1 proB vs pre-pro B
db.ko1 = dba.contrast(db.ko, group1 = db.ko$masks$`Pro-B`,
                      group2 = db.ko$masks$`Pre-pro-B`) %>%
  dba.analyze(., method = DBA_ALL_METHODS)
## ko2 small pre B vs proB
db.ko2 = dba.contrast(db.ko, group1 = db.ko$masks$`small pre-B` ,
                      group2 = db.ko$masks$`Pro-B`) %>%
  dba.analyze(., method = DBA_ALL_METHODS)

dba.show(db.ko1, bContrasts = T)
dba.show(db.ko2, bContrasts = T)

dar.ko1.report=dba.report(db.ko1)

write.csv(data.frame(dar.ko1.report),"output/dar.ko1.csv")

export.bed(dar.ko1.report, 'output/dar.ko1.bed')
export.bed(dar.ko1.report %>% filter(Fold>0), 'output/dar.ko1_up.bed')
export.bed(dar.ko1.report %>% filter(Fold<0), 'output/dar.ko1_down.bed')

saveRDS(dar.ko1.report,"output/dar.ko1.report.rds")

dar.ko1.anno <- comp_meth_anno(dar.ko1.report, cpgs_info = NULL, column = "Fold")
ggsave("output/dar.ko1_annotation_pie.pdf", width = 10, height = 7)

dar.ko2.report=dba.report(db.ko2)

write.csv(data.frame(dar.ko2.report),"output/dar.ko2.csv")

export.bed(dar.ko2.report, 'output/dar.ko2.bed')
export.bed(dar.ko2.report %>% filter(Fold>0), 'output/dar.ko2_up.bed')
export.bed(dar.ko2.report %>% filter(Fold<0), 'output/dar.ko2_down.bed')

saveRDS(dar.ko2.report,"output/dar.ko2.report.rds")

dar.ko2.anno <- comp_meth_anno(dar.ko2.report, cpgs_info = NULL, column = "Fold")
ggsave("output/dar.ko2_annotation_pie.pdf", width = 10, height = 7)

## load all peaks as gr
gr_alllist <- lapply(read.csv("samplesheet_atac_ZD.csv")$Peaks %>% as.character, 
                     function(x){
                                  import(x, 
                                         format = "BED",
                                         extraCols = c(singnalValue = "numeric", pValue = "numeric",
                                                       qValue = "numeric", peak = "integer")) %>% reduce_ranges()
                                })

## peak width and number of samples
gr_list <- list(ppb.wt=gr_alllist[1:4], ppb.ko=gr_alllist[5:8], 
                pb.wt=gr_alllist[9:12], pb.ko=gr_alllist[13:16],
                spb.wt=gr_alllist[17:20], spb.ko=gr_alllist[21:24])


cpgs_info <- build_annotations(genome = "mm10", annotations = "mm10_cpgs")

## number of peaks
peak_nums <- sapply(gr_list, lengths, simplify = T,USE.NAMES = T)
## mean peak width
peak_widths <- sapply(gr_list, function(x) {
  sapply(x, function(y) width(y) %>% mean)
})
## mean peak width in CpG islands
peak_cgiwidths <- sapply(gr_list, function(x) {
  sapply(x, function(y) filter_by_overlaps(y, cpgs_info %>% filter(type=="mm10_cpg_islands")) %>% width() %>% mean)
})
## mean peak width in CpG open seas
peak_seawidths <- sapply(gr_list, function(x) {
  sapply(x, function(y) filter_by_overlaps(y, cpgs_info %>% filter(type=="mm10_cpg_inter")) %>% width() %>% mean)
})
## total opened bases
peak_cov <- sapply(gr_list, function(x) {
  sapply(x, function(y) width(y) %>% sum)
})
## opened bases in CpG islands
peak_cgicov <- sapply(gr_list, function(x) {
  sapply(x, function(y) filter_by_overlaps(y, cpgs_info %>% filter(type=="mm10_cpg_islands")) %>% width() %>% sum)
})
## opened bases in CpG open seas
peak_seacov <- sapply(gr_list, function(x) {
  sapply(x, function(y) filter_by_overlaps(y, cpgs_info %>% filter(type=="mm10_cpg_inter")) %>% width() %>% sum)
})

rbind(peak_nums, 
      peak_widths, peak_cgiwidths, peak_seawidths, 
      peak_cov, peak_cgicov, peak_seacov) %>% data.frame() %>%
  mutate(type=rep(c("number","width","cgiwidth","seawidth","coverage","cgicoverage","seacoverage"), each=4) %>% 
           factor(levels = c("width","cgiwidth","seawidth","coverage","cgicoverage","seacoverage","number"))) %>%
  pivot_longer(1:6) %>% 
  separate(name, c("celltype","genotype"), remove=F) %>%
  ggplot(aes(x=name, y=value, fill=celltype)) + 
  geom_boxplot()+
  facet_wrap(.~type, scales = "free_y")+
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("ppb.wt","ppb.ko"),c("pb.wt", "pb.ko"), c("spb.wt","spb.ko")))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  #scale_fill_manual(values = c(hue_pal()(4)[c(3,2,4,1)],"#808080"))
ggsave("output/peak_width_bases_all.pdf", width = 10, height = 10)


## sankey plot

source("join_overlap_full.R")
z <- join_overlap_full(dar.ppb.report %>% mutate(status=ifelse(Fold>0,"open","closed")) %>% plyranges::select(status), 
                       dar.pb.report %>% mutate(status=ifelse(Fold>0,"open","closed") ) %>% plyranges::select(status)) %>% 
  join_overlap_full(., dar.spb.report %>% mutate(status=ifelse(Fold>0,"open","closed") ) %>% plyranges::select(status))
colnames(mcols(z)) <- c("spb","pb","ppb")

z %>% mcols() %>%
  data.frame() %>% 
  # mutate_all(function(fac) {
  #   factor(fac, levels = levels(addNA(fac)), labels = c("down","up","non"), 
  #          exclude = NULL)
  # }) %>%
  mutate_all(function(fac) {
    factor(fac, levels = c("open","closed", NA), labels = c("up","down","non"), 
           exclude = NULL)
  }) %>%
  count(ppb, pb, spb) %>%
  to_lodes_form(axes=1:3, id="id") %>%
  ggplot(aes(x=x, y=n, stratum=stratum, alluvium=id,
             fill=stratum, label=stratum)) +
  geom_flow(stat = "alluvium", color="darkgray")+
  geom_stratum(alpha=0.5) + 
  scale_fill_manual(values = c("#00BA38","#F8766D","#619CFF")) +
  ylab("Number of regions") + xlab("KO vs WT")
ggsave("output/atac_sankey_all.pdf", width = 4, height = 4)

## pathway analysis
source("~/ls25_scratch/feng.yan/SA_ATAC_chemo/atac_diff_SA_all/code/pathways_analysis.R")
dar.pb.path <- pathways_analysis(dar.pb.report, dir = "output/dar.pb")
dar.ppb.path <- pathways_analysis(dar.ppb.report, dir = "output/dar.ppb")
dar.spb.path <- pathways_analysis(dar.spb.report, dir = "output/dar.spb")
dar.wt1.path <- pathways_analysis(dar.wt1.report, dir = "output/dar.wt1")
dar.wt2.path <- pathways_analysis(dar.wt2.report, dir = "output/dar.wt2")
dar.ko1.path <- pathways_analysis(dar.ko1.report, dir = "output/dar.ko1")
dar.ko2.path <- pathways_analysis(dar.ko2.report, dir = "output/dar.ko2")

## heatmap of top 5k in each KO vs WT

merged.dar <- join_overlap_full(dar.ppb.report %>% head(n=5000),  # .x
                                dar.pb.report %>% head(n=5000)) %>%  # .y
  join_overlap_full(., dar.spb.report %>% head(n=5000)) # 
mcols(merged.dar) <- NULL

idx <- findOverlaps(merged.dar, peaks) %>% subjectHits()

anno.col <- data.frame(condition = db$samples$Condition %>% factor(., levels = c('WT','Ki67KO')), 
                       tissue = db$samples$Tissue,
                       row.names = db$samples$SampleID)

anno.row <- data.frame(ppb = ifelse(merged.dar %in% dar.ppb.report,'DAR','Non'),
                       pb = ifelse(merged.dar %in% dar.pb.report,'DAR','Non'),
                       spb = ifelse(merged.dar %in% dar.spb.report,'DAR','Non'),
                      row.names = idx)

anno.color <- list(tissue = color.celltype,
                   condition = hue_pal()(2)[c(2,1)],
                   pb = c(color.celltype[2], '#56B4E9'),
                   ppb = c(color.celltype[1], '#56B4E9'),
                   spb = c(color.celltype[3], '#56B4E9'))
names(anno.color$tissue) <- c('Pre-pro-B','Pro-B','small pre-B')
names(anno.color$condition) <- c('WT','Ki67KO')
names(anno.color$pb) <- c('DAR','Non')
names(anno.color$ppb) <- c('DAR','Non')
names(anno.color$spb) <- c('DAR','Non')
  
db_mat <- db$binding[idx, 4:27] %>% log2() 

pheatmap(db_mat, scale = 'row',
           cluster_rows = F,
           annotation_row = anno.row, 
           annotation_col = anno.col,
           annotation_colors = anno.color,
           show_rownames = F,
         clustering_callback = function(hc, mat){
           sv = svd(t(mat))$v[,1]
           dend = reorder(as.dendrogram(hc), wts = rep(c(1,1,1,1,2,2,2,2), 3))
           as.hclust(dend)
         },
           filename = 'plot/heatmap_all_top5k.pdf',
           height = 6, width = 5
           )


## summarize annotation by CpG island and tss

anno2 <- peaks %>% comp_meth_anno(., cpgs_info = cpgs_info, simplify = T)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
anno2.all <- rbind(dar.pb.anno$summary %>% data.frame %>% mutate(comparison = 'pb'), 
                   dar.ppb.anno$summary %>% data.frame %>% mutate(comparison = 'ppb'), 
                   dar.spb.anno$summary %>% data.frame %>% mutate(comparison = 'spb'),
                   dar.wt1.anno$summary %>% data.frame %>% mutate(comparison = 'pb2ppb'), 
                   dar.wt2.anno$summary %>% data.frame %>% mutate(comparison = 'pb2spb')) %>%
  data.frame() %>%
  rownames_to_column() %>% 
  mutate(category = ifelse(str_detect(rowname, 'mm10'), 'cpg','tss'),
         location = gsub('[0-9]$','',rowname)) %>% dplyr::select(-rowname)

lapply(c("fill","stack"), function(x){
  y <- ifelse(x=="fill","pct","num")
  ylab <- ifelse(x=="fill", "Proportion (%)", "Numbers of DARs")
  p <- anno2.all %>% dplyr::select(c(2,3,7:9)) %>% dplyr::filter(category=="cpg") %>%
    gather(key = "group", value = "number", 1:2) %>% 
    mutate(CGIanno=factor(location, 
                          levels = paste("mm10", c("cpg_inter","cpg_shelves","cpg_shores", "cpg_islands"), sep = "_"),
                          labels = c("OpenSea","Shelves","Shores", "Islands"))) %>%
    tidyr::unite("x", comparison, group, remove = F) %>% 
    # mutate(x=factor(x, levels = c("wl2lwe_up", "wl2lwe_down", "le2we_up","le2we_down","ll2le_up", "ll2le_down",
    #                               "tall2ll_up", "tall2ll_down"))) %>% 
    ggplot(aes(x = fct_rev(x), y = number, fill = CGIanno)) + 
    geom_bar(stat = "identity", position = x) + 
    scale_fill_manual(values=cbPalette[c(1,2,8,7)]) +
    coord_flip() + 
    ylab(ylab) + xlab("Comparison") 
  ggsave(paste0("plot/anno_all_bar_",y,".pdf"), p, width = 5, height = 3)
})

lapply(c("fill","stack"), function(x){
  y <- ifelse(x=="fill","pct","num")
  ylab <- ifelse(x=="fill", "Proportion (%)", "Numbers of DARs")
  p <- anno2$summary %>% data.frame() %>%
    mutate(comparison=rep("allcpg", each=13),
           category=c(rep("cpg",4),(rep("tss",9))),
           location=rownames(anno2$summary)) %>%
    filter(category=="cpg") %>%
    dplyr::select(-2) %>% 
    mutate(CGIanno=factor(location, 
                          levels = paste("mm10", c("cpg_inter","cpg_shelves","cpg_shores", "cpg_islands"), sep = "_"),
                          labels = c("OpenSea","Shelves","Shores", "Islands"))) %>%
    ggplot(aes(x = comparison, y = total, fill = CGIanno)) + 
    geom_bar(stat = "identity", position = x) + 
    scale_fill_manual(values=cbPalette[c(1,2,8,7)]) +
    coord_flip() + 
    ylab(ylab) + xlab("background") 
  ggsave(paste0("plot/anno_background_bar_",y,".pdf"), p, width = 5, height = 1)
})

lapply(c("fill","stack"), function(x){
  y <- ifelse(x=="fill","pct","num")
  ylab <- ifelse(x=="fill", "Proportion (%)", "Numbers of DARs")
  p <- anno2.all %>% dplyr::select(c(2,3,7:9)) %>% dplyr::filter(category=="tss") %>%
    gather(key = "group", value = "number", 1:2) %>% 
    mutate(tssanno=factor(location, 
                          levels = anno2.all$location[5:13],
                          labels = c("Other","Other","Distal", "Other",'Exon','Intron','Promoter','Promoter','Promoter'))) %>%
    tidyr::unite("x", comparison, group, remove = F) %>% 
    ggplot(aes(x = fct_rev(x), y = number, fill = tssanno)) + 
    geom_bar(stat = "identity", position = x) + 
    scale_fill_manual(values=cbPalette[c(3,1,2,8,7)]) +
    coord_flip() + 
    ylab(ylab) + xlab("Comparison") 
  ggsave(paste0("plot/anno_all_bar_tss_",y,".pdf"), p, width = 5, height = 3)
})

lapply(c("fill","stack"), function(x){
  y <- ifelse(x=="fill","pct","num")
  ylab <- ifelse(x=="fill", "Proportion (%)", "Numbers of DARs")
  p <- anno2$summary %>% data.frame() %>%
    mutate(comparison=rep("allcpg", each=13),
           category=c(rep("cpg",4),(rep("tss",9))),
           location=rownames(anno2$summary)) %>%
    filter(category=="tss") %>%
    dplyr::select(-2) %>% 
    mutate(tssanno=factor(location, 
                          levels = anno2.all$location[5:13],
                          labels = c("Other","Other","Distal", "Other",'Exon','Intron','Promoter','Promoter','Promoter'))) %>%
    ggplot(aes(x = comparison, y = total, fill = tssanno)) + 
    geom_bar(stat = "identity", position = x) + 
    scale_fill_manual(values=cbPalette[c(3,1,2,8,7)]) +
    coord_flip() + 
    ylab(ylab) + xlab("background") 
  ggsave(paste0("plot/anno_background_bar_tss_",y,".pdf"), p, width = 5, height = 1)
})


dar.seq <- lapply(list(dar.pb.report, dar.ppb.report, dar.spb.report,
            dar.wt1.report, dar.wt2.report), function(x){
              cts <- count(x %>% data.frame(), seqnames)
              width <- x %>% data.frame() %>% group_by(seqnames) %>%
                summarise(width = sum(width))
              full_join(cts, width, by='seqnames')
            }) %>% reduce(full_join, by='seqnames') %>%
  mutate(seqnames = factor(seqnames, levels = paste0('chr',c(1:19,'X','Y'))))
colnames(dar.seq) <- c('chr',paste(rep(c('pb','ppb','spb','pb2ppb','pb2spb'), each=2), 
                          rep(c('num','width'), 5), sep = '_'))

dar_num <- dar.seq %>% pivot_longer(!chr) %>%
  separate(name, c('comparison', 'feature')) %>% 
  filter(feature == 'num' & (comparison %in% c('pb','ppb','spb'))) 
write.csv(dar_num, 'output/number_dars.csv')

dar_num %>% 
  ggplot(aes(x=chr, y=value, group=comparison, col=comparison)) + 
  geom_point() + geom_line() +
  ylab('Number of DARs') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('plot/number_dars_by_chr.pdf', width = 10, height = 5)

seqlen <- getChromInfoFromUCSC('mm10') %>% 
  dplyr::filter(chrom %in% levels(dar.seq$chr))

dar_pct <- dar.seq %>% left_join(seqlen, by=c('chr'='chrom')) %>% 
  left_join(peaks %>% data.frame() %>% group_by(seqnames) %>%
              summarise(width = sum(width)), by=c('chr'='seqnames')) %>%
  pivot_longer(pb_num:pb2spb_width, values_to = 'dar_bases') %>%
  separate(name, c('comparison', 'feature')) %>% 
  filter(feature == 'width' & (comparison %in% c('pb','ppb','spb'))) %>% 
  rename(chromsize = length, atac_background = width) %>%
  mutate(pct_chr = dar_bases / chromsize * 100, 
         pct_atac = dar_bases / atac_background *100) %>% 
  pivot_longer(starts_with('pct'), values_to = 'percentage')
write.csv(dar_pct, 'output/percentage_covered_dars.csv')

dar_pct %>%
  ggplot(aes(x=chr, y=percentage, group=comparison, col=comparison)) + 
  geom_point() + geom_line() +
  ylab('% altered') +
  facet_wrap(.~name, scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('plot/cov_dars_by_chr.pdf', width = 10, height = 5)


write.csv(data.frame(dar.ppb.anno$tss),"output/dar.ppb_anno.csv")
write.csv(data.frame(dar.pb.anno$tss),"output/dar.pb_anno.csv")
write.csv(data.frame(dar.spb.anno$tss),"output/dar.spb_anno.csv")