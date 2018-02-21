## figures
library(ggplot2)
library(DEGreport)
load("data/malog.rda")
keep = rownames(meta)[!grepl("re", rownames(meta), ignore.case = TRUE)]
meta = meta[keep,]
malog = malog[,keep]

#note: always pass alpha on the 0-255 scale
makeTransparent<-function(someColor, alpha=100)
{
    newColor<-col2rgb(someColor)
    apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                                blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

cc = RColorBrewer::brewer.pal(n = 8, name = "RdYlBu")
cc[2:7] = makeTransparent(cc[2:7], 50)

degPCA(malog, meta, condition = "Timepoint") + 
    scale_color_manual(values = cc) +
    ggsave("figures/pca_pc1_pc2.pdf")
degPCA(malog, meta, condition = "Timepoint", pc1 = "PC1", pc2 = "PC4") +
    scale_color_manual(values = cc) +
    ggsave("figures/pca_pc1_pc4.pdf")

keep2 = rownames(meta)[grepl("Fast_D1", rownames(meta))]
degPCA(malog[,keep2], meta[keep2,], condition = "Timepoint", name = "Description")


## covariates
load("data/malog.rda")
library(readr)
library(dplyr)
library(limma)
profile = read_tsv("clean_timepoint_de.xls") %>% 
    filter(!is.na(cluster))

c = c("NEFA","Ketones","Insulin","Glucagon","Leptin",
  "FGF21","T3","Adiponectin","HOMA-IR")
h = c("Cortisol","21-Deoxycortisol","Cortexolone","Cortisone")

m = meta[,c]
h = malog[h,rownames(m)]

new_meta = bind_cols(m, as.data.frame(t(h)))
rownames(new_meta) = rownames(m)

model <- model.matrix(data=meta, ~ Subject + Timepoint)
v = vooma(malog, design = model, plot = T)

cor_meta <- lapply(unique(profile$cluster), function(c){
    .genes <- as.character(profile$gene[profile$cluster==c])
    meds = colMedians(v$E[.genes,])
    lapply(colnames(new_meta), function(m){
        .meta_noempty = as.numeric(as.factor(new_meta[!is.na(new_meta[,m]),m]))
        .cor = cor.test(.meta_noempty, meds[!is.na(new_meta[,m])])
        data.frame(covar=m, pval=.cor$p.value, rho=.cor$estimate, cluster=c)
    }) %>% bind_rows()
    
}) %>% bind_rows()

cor_meta %>% 
    mutate(fdr=p.adjust(pval, method="fdr")) %>% 
    write_tsv("clean_cor_meta_13vars.xls")
