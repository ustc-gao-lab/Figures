library(ComplexHeatmap)
library(circlize)
library(clipr)
library(matrixStats)

setwd("C:/Users/gaoda/OneDrive - USTC/Gao lab/Scripts/metabolic/")
sample = read.table('sample.txt', sep='\t', header=FALSE, row.names=1)
colnames(sample) = c('genotype', 'sex')
sample_df = data.frame(genotype=sample$genotype, condition=sample$sex)
sample_annot = HeatmapAnnotation(df=sample_df, gap=unit(1,"mm"), width=unit(5,"mm"), 
                                 col=list(genotype=c('WT'='#00cc99','LCK-cKI'='#66ccff'),
                                          condition=c('male'='#4d4d4d','female'='#e6e6e6')))
MS_result=read_clip_tbl(sep
                        ="\t")
#MS_result=MS_result[2:nrow(MS_result),]
rownames_MS_result=MS_result[[1]]
MS_result[,c(1)]=NULL
colnames(MS_result)=c("WT_male_rep1","WT_male_rep2","LCK_cKI_male_rep1","LCK_cKI_male_rep2","WT_female_rep1","WT_female_rep2","LCK_cKI_female_rep1","LCK_cKI_female_rep2")
MS_result=as.data.frame(sapply(as.data.frame(MS_result),as.numeric))
rownames(MS_result)=rownames_MS_result
MS_result[is.na(MS_result)]=0
MS_result_Z=MS_result
write.table(MS_result, "MS_result.txt", sep="\t", quote=FALSE)
for(i in 1:8){
  MS_result_Z[[i]]=(MS_result[[i]]-rowMeans(MS_result[,1:8]))/rowSds(as.matrix(MS_result[,1:8]))
}


MS_result_Z[is.na(MS_result_Z)]=0
MS_result_Z[MS_result_Z==Inf]=0

data_heatmap = Heatmap(MS_result_Z, name='z-score', col=colorRamp2(c(-3,0,3),c('#0000ff','#ffffff','#ff6600')),
                       row_title='metabolite changes', row_title_side='left', column_title_side='bottom',column_order=1:8,
                       clustering_method_rows='complete', clustering_distance_rows='manhattan', row_dend_width=unit(10,'mm'),
                       clustering_method_columns='complete', clustering_distance_columns='euclidean', column_dend_height=unit(20,'mm'),
                       show_row_names=FALSE, show_column_names=TRUE, column_names_gp=gpar(fontsize=9), show_heatmap_legend=TRUE, top_annotation=sample_annot)
png(file = "C:/Users/gaoda/OneDrive - USTC/Gao lab/Scripts/metabolic/test_New4.png", width=9, height=9, units='in', res=300)
draw(data_heatmap)
dev.off()


