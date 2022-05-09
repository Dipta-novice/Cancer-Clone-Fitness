#BiocManager::install("infercnv",version = "1.8")
library(infercnv)

# controlData = as.data.frame(t(data.table::fread("C:\Users\acer\OneDrive\Documents\infercnv\GBM-009-01-1C.gene.expression.matrix.tsv\GBM-009-01-1C.gene.expression.matrix.tsv",header = TRUE)))
# colnames(controlData) = controlData[1,]
# controlData = controlData[-1,]
# saveRDS("C:\Users\acer\OneDrive\Documents\infercnv\ControlData5000.rds")

# load("C:\Users\acer\OneDrive\Documents\infercnv\ControlData5000.RData")
# 
# controlData5000 = controlData[,sample(1:ncol(controlData), 5000) ]
# saveRDS(C:\Users\acer\OneDrive\Documents\infercnv\ControlData5000.rds")

controlData = readRDS("C:\Users\acer\OneDrive\Documents\infercnv\ControlData5000.rds")

QueryData = data.matrix(read.table("C:\Users\acer\OneDrive\Documents\infercnv\GBM-009-01-1C.gene.expression.matrix.tsv\GBM-009-01-1C.gene.expression.matrix.tsv"))

Commongenes = intersect(rownames(controlData),rownames(QueryData))
combinedData = cbind(controlData[Commongenes,],QueryData[Commongenes,])

annotations_file = data.frame(cbind(c(rep("Control",length(colnames(controlData))),
                                      rep("Query",length(colnames(QueryData)))),
                                    c(colnames(controlData),colnames(QueryData))))
rownames(annotations_file) = annotations_file$X2
annotations_file = annotations_file[,"X1" ,drop=FALSE]
colnames(annotations_file) = NULL
load("~/Clonal_fitness/gencode_v19_gene_pos.RData")

# Creation of an infercnv object
infercnv_obj1.1 = CreateInfercnvObject(raw_counts_matrix= data.matrix(combinedData),
                                       annotations_file = annotations_file,
                                       delim="\t",
                                       gene_order_file=gencode_v19_gene_pos,
                                       ref_group_names = c("Control"))


infercnv_obj1.2 = infercnv::run(infercnv_obj1.1,
                                cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                out_dir="C:\Users\acer\OneDrive\Documents\infercnv", 
                                cluster_by_groups=FALSE, 
                                denoise=TRUE,
                                HMM=TRUE)

plot_cnv(infercnv_obj1.2,
         out_dir= "C:\Users\acer\OneDrive\Documents\infercnv",
         obs_title="Query",
         ref_title="Control",
         cluster_by_groups=FALSE,
         x.center=1,
         x.range="auto",
         hclust_method='ward.D',
         color_safe_pal=FALSE,
         output_filename = GBM_query,
         output_format="pdf",
         png_res=300,
         dynamic_resize=0
)