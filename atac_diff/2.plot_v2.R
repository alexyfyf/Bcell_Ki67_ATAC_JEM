
library(amap) ## load Dist function
library(RColorBrewer)

## Fig 5A/B (PCA plots)
color.celltype <- brewer.pal(8, 'Dark2')[c(8,4,7)]

gr.pca <- rowSds(db$binding[,3+(1:24)]) %>% order(decreasing = T) %>% 
  head(n=5000)

autoplot(db$binding[gr.pca, 3+(1:24)] %>% t() %>% log2() %>% 
           prcomp(scale.=T) ,
         data=db$samples[1:24,], label=F) + 
  geom_point(aes(col=Treatment, shape=Condition), size=3) +
  scale_color_manual(values = color.celltype)
ggsave("plot2/pca_atac_all_5k.pdf", width = 4, height = 3, units = "in")

autoplot(db$binding[gr.pca, 3+(1:24)] %>% t() %>% log2() %>% 
           prcomp(scale.=T) ,
         data=db$samples[1:24,], label=F, x = 3, y = 4) + 
  geom_point(aes(col=Treatment, shape=Condition), size=3) +
  scale_color_manual(values = color.celltype)
ggsave("plot2/pca_atac_all_5k_pc34.pdf", width = 4, height = 3, units = "in")

autoplot(db$binding[gr.pca, 3+(1:24)] %>% t() %>% log2() %>% 
           prcomp(scale.=T) ,
         data=db$samples[1:24,], label=F, x = 2, y = 3) + 
  geom_point(aes(col=Treatment, shape=Condition), size=3) +
  scale_color_manual(values = color.celltype)
ggsave("plot2/pca_atac_all_5k_pc23.pdf", width = 4, height = 3, units = "in")

## use all peaks
autoplot(db$binding[, 3+(1:24)] %>% t() %>% log2() %>% 
           prcomp(scale.=T) ,
         data=db$samples[1:24,], label=F) + 
  geom_point(aes(col=Treatment, shape=Condition), size=3) +
  scale_color_manual(values = color.celltype)
ggsave("plot2/pca_atac_all_allpeak.pdf", width = 4, height = 3, units = "in")

autoplot(db$binding[, 3+(1:24)] %>% t() %>% log2() %>% 
           prcomp(scale.=T) ,
         data=db$samples[1:24,], label=F, x = 3, y = 4) + 
  geom_point(aes(col=Treatment, shape=Condition), size=3) +
  scale_color_manual(values = color.celltype)
ggsave("plot2/pca_atac_all_allpeak_pc34.pdf", width = 4, height = 3, units = "in")

autoplot(db$binding[, 3+(1:24)] %>% t() %>% log2() %>% 
           prcomp(scale.=T) ,
         data=db$samples[1:24,], label=F, x = 2, y = 3) + 
  geom_point(aes(col=Treatment, shape=Condition), size=3) +
  scale_color_manual(values = color.celltype)
ggsave("plot2/pca_atac_all_allpeak_pc23.pdf", width = 4, height = 3, units = "in")

## separate cell types
autoplot(db$binding[gr.pca, 3+(1:8)] %>% t() %>% log2() %>% 
           prcomp(scale.=T) ,
         data=db$samples[1:8,], label=F) + 
  geom_point(aes(col=Treatment, shape=Condition), size=3) +
  scale_color_manual(values = color.celltype[1])
ggsave("plot2/pca_atac_ppb_5k.pdf", width = 4, height = 3, units = "in")

autoplot(db$binding[gr.pca, 3+(9:16)] %>% t() %>% log2() %>% 
           prcomp(scale.=T) ,
         data=db$samples[9:16,], label=F) + 
  geom_point(aes(col=Treatment, shape=Condition), size=3) +
  scale_color_manual(values = color.celltype[2])
ggsave("plot2/pca_atac_pb_5k.pdf", width = 4, height = 3, units = "in")

autoplot(db$binding[gr.pca, 3+(17:24)] %>% t() %>% log2() %>% 
           prcomp(scale.=T) ,
         data=db$samples[17:24,] , label=F) + 
  geom_point(aes(col=Treatment, shape=Condition), size=3) +
  scale_color_manual(values = color.celltype[3])
ggsave("plot2/pca_atac_spb_5k.pdf", width = 4, height = 3, units = "in")


## Fig 5C (heatmap redo)

save_pheatmap_png <- function(x, filename, width=1000, height=1200, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

dar.all <- full_join(dar.ppb.report %>% data.frame %>% 
  rename(Fold_ppb=Fold, FDR_ppb=FDR) %>%
  unite(region, 1:3, remove = F) %>% dplyr::select(c(1:4,6,10,12)) ,
  dar.pb.report %>% data.frame %>% 
    rename(Fold_pb=Fold, FDR_pb=FDR) %>%
    unite(region, 1:3, remove = F) %>% dplyr::select(c(1:4,6,10,12)) ) %>%
  full_join(dar.spb.report %>% data.frame %>% 
              rename(Fold_spb=Fold, FDR_spb=FDR) %>%
              unite(region, 1:3, remove = F) %>% dplyr::select(c(1:4,6,10,12))) %>% 
  full_join(dar.wt1.report %>% data.frame %>% 
              rename(Fold_pb2ppb_wt=Fold, FDR_pb2ppb_wt=FDR) %>%
              unite(region, 1:3, remove = F) %>% dplyr::select(c(1:4,6,10,12))) %>% 
  full_join(dar.wt2.report %>% data.frame %>% 
              rename(Fold_spb2pb_wt=Fold, FDR_spb2pb_wt=FDR) %>%
              unite(region, 1:3, remove = F) %>% dplyr::select(c(1:4,6,10,12))) %>% 
  full_join(dar.ko1.report %>% data.frame %>% 
              rename(Fold_pb2ppb_ko=Fold, FDR_pb2ppb_ko=FDR) %>%
              unite(region, 1:3, remove = F) %>% dplyr::select(c(1:4,6,10,12))) %>% 
  full_join(dar.ko2.report %>% data.frame %>% 
              rename(Fold_spb2pb_ko=Fold, FDR_spb2pb_ko=FDR) %>%
              unite(region, 1:3, remove = F) %>% dplyr::select(c(1:4,6,10,12))) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

dat.all.anno <- dar.all %>% mutate(heatmap = dar.all %in% merged.dar) %>%
  join_overlap_inner(anno2$tss, minoverlap = 2) 

dat.all.anno %>% data.frame() %>% write.csv('plot2/dar.all.anno.csv')

go <- dat.all.anno %>% data.frame() %>% 
  dplyr::filter(Fold_pb <0  & Fold_spb <0) %>% 
  dplyr::pull(geneId) %>% 
  unique() %>% 
  enrichGO(OrgDb = org.Mm.eg.db, ont = 'BP', readable = T) %>%
  simplify() 

go %>%
  cnetplot(showCategory = 5 )

## Fig 5D (barplot to replace sankey, and volcano plots)
dar.barplot <- lapply(list(dar.ppb.report, dar.pb.report, dar.spb.report), function(x){
                         cts <- x %>% data.frame() %>% count(Fold>0)
                         }) %>% reduce(full_join, by='Fold > 0') 
colnames(dar.barplot) <- c('up','ppb','pb','spb')

dar.barplot %>% pivot_longer(!up) %>%
  mutate(name=factor(name, levels=c('ppb','pb','spb'))) %>%
  mutate(DAR=ifelse(up,'Up','Down')) %>%
  ggplot(aes(x=name, y=value, fill=DAR)) + 
  geom_bar(stat='identity') +
  ylab('Number of DARs') +
  scale_fill_manual(values = c("#F8766D","#00BA38"))
ggsave('plot2/number_dar.pdf', width = 3, height = 2.5)


df.ppb <- data.frame(dar.ppb.anno$tss) %>%
  mutate(DAR = ifelse(FDR>1e-5, 'A',ifelse(Fold>1, 'Up',ifelse(Fold< -1, 'Down', 'B'))),
         log10FDR = -log10(FDR))
df.ppb %>%
  ggplot(aes(x=Fold, y=log10FDR, col=DAR)) + 
  geom_point() + 
  geom_hline(yintercept = 5, linetype='dashed', col = 'grey') + 
  geom_vline(xintercept = c(-1,1), linetype='dashed', col = 'grey') +
  # geom_text_repel(data = df.ppb %>%
  #                   filter(FDR<1e-5 & abs(Fold) > 1), # labels only genes in the interesting pathway
  #                 aes(label = SYMBOL ), #, color = sig),
  #                 size = 3.5,
  #                 color = "black",
  #                 nudge_x = 0.3, nudge_y = 0.1, max.overlaps = 50) +
  scale_color_manual(values = c('#666666','#666666' , "#F8766D","#00BA38")) +
  ggtitle('Volcano Plot for Pre-Pro B KO vs WT (Up 158, Down 320)')
ggsave('plot2/vp_ppb.pdf', width = 6, height = 6)


df.pb <- data.frame(dar.pb.anno$tss) %>%
  mutate(DAR = ifelse(FDR>1e-5, 'A',ifelse(Fold>1, 'Up',ifelse(Fold< -1, 'Down', 'B'))),
         log10FDR = -log10(FDR))
df.pb %>%
  ggplot(aes(x=Fold, y=log10FDR, col=DAR)) + 
  geom_point() + 
  geom_hline(yintercept = 5, linetype='dashed', col = 'grey') + 
  geom_vline(xintercept = c(-1,1), linetype='dashed', col = 'grey') +
  # geom_text_repel(data = df.pb %>%
  #                   filter(FDR<1e-5 & abs(Fold) > 1), # labels only genes in the interesting pathway
  #                 aes(label = SYMBOL ), #, color = sig),
  #                 size = 3.5,
  #                 color = "black",
  #                 nudge_x = 0.3, nudge_y = 0.1, max.overlaps = 50) +
  scale_color_manual(values = c('#666666','#666666' , "#F8766D","#00BA38")) +
  ggtitle('Volcano Plot for Pro B KO vs WT (Up 3588, Down 8436)')
ggsave('plot2/vp_pb.pdf', width = 6, height = 6)

df.spb <- data.frame(dar.spb.anno$tss) %>%
  mutate(DAR = ifelse(FDR>1e-5, 'A',ifelse(Fold>1, 'Up',ifelse(Fold< -1, 'Down', 'B'))),
         log10FDR = -log10(FDR))
df.spb %>%
  ggplot(aes(x=Fold, y=log10FDR, col=DAR)) + 
  geom_point() + 
  geom_hline(yintercept = 5, linetype='dashed', col = 'grey') + 
  geom_vline(xintercept = c(-1,1), linetype='dashed', col = 'grey') +
  # geom_text_repel(data = df.spb %>%
  #                   filter(FDR<1e-5 & abs(Fold) > 1), # labels only genes in the interesting pathway
  #                 aes(label = SYMBOL ), #, color = sig),
  #                 size = 3.5,
  #                 color = "black",
  #                 nudge_x = 0.3, nudge_y = 0.1, max.overlaps = 50) +
  scale_color_manual(values = c('#666666','#666666' , "#F8766D","#00BA38")) +
  ggtitle('Volcano Plot for Small-Pre B KO vs WT (Up 6996, Down 2289)')
ggsave('plot2/vp_spb.pdf', width = 6, height = 6)


## Fig 5E (percentage of chromosome altered redo)
dar_pct %>%
  filter(name=='pct_chr', comparison != 'ppb') %>%
  mutate(chr = factor(chr, levels = paste0('chr', c(1:19, 'X','Y')))) %>%
  ggplot(aes(x=chr, y=percentage, group=comparison, col=comparison)) + 
  geom_point() + geom_line() +
  ylab('% altered') +
  facet_wrap(.~name, scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values = color.celltype[2:3])
  
ggsave('plot2/cov_dars_by_chr.pdf', width = 10, height = 5)

library(treemapify)
library(viridis)

dar_pct %>%
  filter(name=='pct_chr', comparison == 'pb') %>%
  filter(!(chr %in% c('chrY','chrX'))) %>%
  ggplot(aes(area = percentage, fill = percentage, label = chr )) +
  geom_treemap() +
  geom_treemap_text(colour = "white",
                    place = "centre",
                    size = 15) +
  scale_fill_viridis_c()
ggsave('plot2/pb_pct_chr_mosaic.pdf', width = 5, height = 4)

dar_pct %>%
  filter(name=='pct_chr', comparison == 'spb') %>%
  filter(!(chr %in% c('chrY','chrX'))) %>%
  ggplot(aes(area = percentage, fill = percentage, label = chr )) +
  geom_treemap() +
  geom_treemap_text(colour = "white",
                    place = "centre",
                    size = 15) +
  scale_fill_viridis_c()
ggsave('plot2/spb_pct_chr_mosaic.pdf', width = 5, height = 4)

dar_pct %>%
  filter(name=='pct_atac', comparison == 'pb') %>%
  filter(!(chr %in% c('chrY','chrX'))) %>%
  ggplot(aes(area = percentage, #fill = percentage, 
             label = paste(chr, paste0(round(percentage,2), '%'), sep = '\n') ) ) +
  geom_treemap() +
  geom_treemap_text(colour = "white",
                    place = "centre",
                    size = 15) +
  scale_fill_viridis_c()
ggsave('plot2/pb_pct_atac_mosaic.pdf', width = 5, height = 4)

dar_pct %>%
  filter(name=='pct_atac', comparison == 'spb') %>%
  filter(!(chr %in% c('chrY','chrX'))) %>%
  ggplot(aes(area = percentage, #fill = percentage, 
             label = paste(chr, paste0(round(percentage,2), '%'), sep = '\n')  )) +
  geom_treemap() +
  geom_treemap_text(colour = "white",
                    place = "centre",
                    size = 15) +
  scale_fill_viridis_c()
ggsave('plot2/spb_pct_atac_mosaic.pdf', width = 5, height = 4)

## sankey compare DARs in WT and Ki67
source('~/github/utils/plot_sankey.R')
plot_sankey(list(WT=dar.wt1.report, KO=dar.ko1.report), column = 'Fold', 
            fulljoin = T)$plot + 
  ylab("Number of regions") + xlab("DARs ProB vs PreProB")
ggsave("plot2/atac_sankey_pb2ppb.pdf", width = 4, height = 4)

plot_sankey(list( WT=dar.wt2.report, KO=dar.ko2.report), column = 'Fold', 
            fulljoin = T)$plot + 
  ylab("Number of regions") + xlab("DARs SmallPreB vs ProB")
ggsave("plot2/atac_sankey_spb2pb.pdf", width = 4, height = 4)

plot_sankey(list(KOP2WTP=dar.pb.report, KOSP2KOP=dar.ko2.report, 
                 WTSP2WTP=dar.wt2.report, KOSP2WTSP=dar.spb.report), 
            column = 'Fold', 
            fulljoin = T) + 
  ylab("Number of regions") + xlab("DARs ProB vs PreProB")

## volcano plot for HOMER results
library(data.table) 
library(EnhancedVolcano)

cols <- c( "down" = "#F8766D", "non" = "#666666", "up" = "#00BA38")
sizes <- c("down" = 2, "up" = 2, "non" = 1) 

homer.ppb <- rbindlist(list(readxl::read_excel('../New_Fig6/homer_summary.xlsx', sheet = 2) %>% mutate(change='down'),
                   readxl::read_excel('../New_Fig6/homer_summary.xlsx', sheet = 3) %>% mutate(change='up')),
                   use.names=FALSE) %>%
  dplyr::select(c(1,4,7,9,10)) %>% 
  dplyr::rename(motif = `Motif Name`, logP = `Log P-value`,
                pct_target = `% of Target Sequences with Motif`, 
                pct_bg = `% of Background Sequences with Motif`) %>%
  filter(pct_target > 0 & pct_bg > 0) %>%
  mutate(enrichment = pct_target/pct_bg) %>%
  filter(enrichment > 1 & !grepl('Unknown', motif)) %>% 
  mutate(log2odd = ifelse(change == 'up', log2(enrichment), -log2(enrichment)),
         sig = ifelse(log2odd > 1 & logP < -5, 'up', 
                      ifelse(log2odd < -1 & logP < -5, 'down', 
                             'non'))) %>%
  separate(motif, c('TF','geo','source'), remove = F, sep = '/') %>%
  separate(TF, c('tf','family'), sep = '\\(|\\)')

ggplot(homer.ppb, aes(x=log2odd, y=-logP)) + 
  geom_point(aes(colour = sig, 
                 size = sig)) + 
  geom_hline(yintercept = c(5), linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  xlab("log2 Odd Ratio") + ylab("-log P value") + 
  scale_color_manual(values = cols) +
  scale_size_manual(values = sizes) +
  geom_text_repel(data = homer.ppb %>% filter(tf %in% c('IRF8','SpiB','PU.1','PU.1:IRF8','RUNX1','RUNX','PU.1-IRF','ERG','E2A',
                                                        'Bach2')),
                  aes(label = tf),
                  size = 3.5,
                  color = "black",
                  nudge_x = 0.3, nudge_y = 0.1, max.overlaps = 50) +
  geom_point(data = homer.ppb %>% filter(tf %in% c('IRF8','SpiB','PU.1','PU.1:IRF8','RUNX1','RUNX','PU.1-IRF','ERG','E2A',
                                                   'Bach2')),
             color = "black",
             size = 1)
ggsave('plot2/evp_ppb_ko2wt_homer.pdf', width = 6, height = 6)

homer.pb <- rbindlist(list(readxl::read_excel('../New_Fig6/homer_summary.xlsx', sheet = 4) %>% mutate(change='down'),
                            readxl::read_excel('../New_Fig6/homer_summary.xlsx', sheet = 5) %>% mutate(change='up')),
                       use.names=FALSE) %>%
  dplyr::select(c(1,4,7,9,10)) %>% 
  dplyr::rename(motif = `Motif Name`, logP = `Log P-value`,
                pct_target = `% of Target Sequences with Motif`, 
                pct_bg = `% of Background Sequences with Motif`) %>%
  filter(pct_target > 0 & pct_bg > 0) %>%
  mutate(enrichment = pct_target/pct_bg) %>%
  filter(enrichment > 1 & !grepl('Unknown', motif)) %>% 
  mutate(log2odd = ifelse(change == 'up', log2(enrichment), -log2(enrichment)),
         sig = ifelse(log2odd > 1 & logP < -5, 'up', 
                      ifelse(log2odd < -1 & logP < -5, 'down', 
                             'non'))) %>%
  separate(motif, c('TF','geo','source'), remove = F, sep = '/') %>%
  separate(TF, c('tf','family'), sep = '\\(|\\)')

ggplot(homer.pb, aes(x=log2odd, y=-logP)) + 
  geom_point(aes(colour = sig, 
                 size = sig)) + 
  geom_hline(yintercept = c(5), linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  xlab("log2 Odd Ratio") + ylab("-log P value") + 
  scale_color_manual(values = cols) +
  scale_size_manual(values = sizes) +
  geom_text_repel(data = homer.pb %>% filter(tf %in% c('IRF8','SpiB','PU.1','PU.1:IRF8','RUNX1','RUNX','PU.1-IRF','ERG','E2A','IRF4', 'EBF', 'PAX5',
                                                       'Bach2','Bcl6')),
                  aes(label = tf),
                  size = 3.5,
                  color = "black",
                  nudge_x = 0.3, nudge_y = 0.1, max.overlaps = 50) +
  geom_point(data = homer.pb %>% filter(tf %in% c('IRF8','SpiB','PU.1','PU.1:IRF8','RUNX1','RUNX','PU.1-IRF','ERG','E2A','IRF4', 'EBF', 'PAX5',
                                                  'Bach2','Bcl6')),
             color = "black",
             size = 1)
ggsave('plot2/evp_pb_ko2wt_homer.pdf', width = 6, height = 6)

homer.spb <- rbindlist(list(readxl::read_excel('../New_Fig6/homer_summary.xlsx', sheet = 6) %>% mutate(change='down'),
                           readxl::read_excel('../New_Fig6/homer_summary.xlsx', sheet = 7) %>% mutate(change='up')),
                      use.names=FALSE) %>%
  dplyr::select(c(1,4,7,9,10)) %>% 
  dplyr::rename(motif = `Motif Name`, logP = `Log P-value`,
                pct_target = `% of Target Sequences with Motif`, 
                pct_bg = `% of Background Sequences with Motif`) %>%
  filter(pct_target > 0 & pct_bg > 0) %>%
  mutate(enrichment = pct_target/pct_bg) %>%
  filter(enrichment > 1 & !grepl('Unknown', motif)) %>% 
  mutate(log2odd = ifelse(change == 'up', log2(enrichment), -log2(enrichment)),
         sig = ifelse(log2odd > 1 & logP < -5, 'up', 
                      ifelse(log2odd < -1 & logP < -5, 'down', 
                             'non'))) %>%
  separate(motif, c('TF','geo','source'), remove = F, sep = '/') %>%
  separate(TF, c('tf','family'), sep = '\\(|\\)')

ggplot(homer.spb, aes(x=log2odd, y=-logP)) + 
  geom_point(aes(colour = sig, 
                 size = sig)) + 
  geom_hline(yintercept = c(5), linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  xlab("log2 Odd Ratio") + ylab("-log P value") + 
  scale_color_manual(values = cols) +
  scale_size_manual(values = sizes) +
  geom_text_repel(data = homer.spb %>% filter(tf %in% c('IRF8','SpiB','PU.1','PU.1:IRF8','RUNX1','PU.1-IRF','ERG','E2A','IRF4', 'EBF', 'PAX5',
                                                        'Bach2','Bcl6','EBF1','RUNX')),
                  aes(label = tf),
                  size = 3.5,
                  color = "black",
                  nudge_x = 0.3, nudge_y = 0.1, max.overlaps = 50) +
  geom_point(data = homer.spb %>% filter(tf %in% c('IRF8','SpiB','PU.1','PU.1:IRF8','RUNX1','PU.1-IRF','ERG','E2A','IRF4', 'EBF', 'PAX5',
                                                   'Bach2','Bcl6','EBF1','RUNX')),
             color = "black",
             size = 1)
ggsave('plot2/evp_spb_ko2wt_homer.pdf', width = 6, height = 6)

## proportion of DARs in genomic locations
anno2.all %>% dplyr::select(c(2,3,7:9)) %>% 
  dplyr::filter(category=="tss" & comparison %in% c('pb','spb')) %>%
  gather(key = "group", value = "number", 1:2) %>%
  data.frame() %>%
  mutate(tssanno=factor(location, 
                        levels = anno2.all$location[5:13],
                        labels = c("Other","Other","Distal", "Other",'Exon','Intron','Promoter','Promoter','Promoter'))) %>%
  tidyr::unite("x", comparison, group, remove = F) %>%
  ggplot(aes(x = x, y = number, fill = tssanno)) + 
  geom_bar(stat = "identity", position = 'fill') + 
  scale_fill_manual(values=cbPalette[c(3,1,2,8,7)]) +
  coord_flip() + 
  ylab('Percentage (%)') + xlab("DARs") 
ggsave('plot2/annobar_tss_pct.pdf', width = 5, height = 2)

anno2.all %>% dplyr::select(c(2,3,7:9)) %>% 
  dplyr::filter(category=="cpg" & comparison %in% c('pb','spb')) %>%
  gather(key = "group", value = "number", 1:2) %>%
  data.frame() %>%
  mutate(CGIanno=factor(location, 
                        levels = paste("mm10", c("cpg_inter","cpg_shelves","cpg_shores", "cpg_islands"), sep = "_"),
                        labels = c("OpenSea","Shelves","Shores", "Islands"))) %>%  tidyr::unite("x", comparison, group, remove = F) %>%
  ggplot(aes(x = x, y = number, fill = CGIanno)) + 
  geom_bar(stat = "identity", position = 'fill') + 
  scale_fill_manual(values=cbPalette[c(3,1,2,8,7)]) +
  coord_flip() + 
  ylab('Percentage (%)') + xlab("DARs") 
ggsave('plot2/annobar_cpg_pct.pdf', width = 5, height = 2)


# Rag1/2 enhancer: chr2:101608312-101608913

data.frame(atac = db$binding[67634, 4:27],
           Genotype = factor(db$samples$Condition, levels = c('WT','Ki67KO')),
           Celltype = db$samples$Tissue) %>%
  ggplot(aes(x = Genotype, y = atac)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(.~Celltype, scales = 'free') +
  stat_compare_means() + 
  ylab('log2 Normalized Read Counts')
ggsave('plot2/rag1_2_enhancer_acc.pdf', width = 8, height = 3)

