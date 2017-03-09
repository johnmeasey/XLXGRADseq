# PCA

Running this script:

```R

library("devtools")
library(gdsfmt)
library(SNPRelate)



vcf.fn <- "XLXVXGXM_merged_sorted.bam.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")

genofile = openfn.gds("test.gds", readonly=FALSE)

add.gdsn(genofile, "sample.annot", samp.annot)

snpgdsSummary("test.gds")

pca <- snpgdsPCA(genofile, num.thread=2)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))


tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

#plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
#text(tab$EV2, tab$EV1,labels=tab$sample.id, cex= 0.4)

library(ggplot2)
#ggplot(...)+...+ theme(axis.text.x = element_text(angle=60, hjust=1))
#devtools::install_github("slowkow/ggrepel")
library(ggrepel)

pdf("PCA_plot_XLXG.pdf",w=8, h=8, version="1.4", bg="transparent")
tab$Species <- c("victorianus", "victorianus", "victorianus", "victorianus", "victorianus", "victorianus", "victorianus", "victorianus", "laevis", "laevis", "gilli", "gilli", "laevis", "gilli", "gilli", "gilli", "gilli", "gilli", "laevis", "gilli", "gilli", "gilli", "gilli", "gilli", "laevis", "laevis", "laevis", "laevis", "laevis", "laevis", "laevis", "laevis", "laevis", "laevis", "laevis", "laevis", "laevis", "muelleri")
tab$samp.color <- c("gray", "gray", "gray", "gray", "gray", "gray", "gray", "gray", "green", "purple", "red", "red", "blue", "red", "red", "red", "red", "red", "blue", "red", "red", "red", "red", "red", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "black")
tab$samp.fieldid <- c("BJE1488_sorted.bam", "BJE1489_sorted.bam", "BJE261_sorted.bam", "BJE263_sorted.bam", "BJE264_sorted.bam", "BJE265_sorted.bam", "BJE266_sorted.bam", "BJE267_sorted.bam", "BJE3545_sorted.bam", "BJE3639_sorted.bam", "XG12_07_sorted.bam", "XG153_sorted.bam", "XG92_sorted.bam", "XGL713_123_sorted.bam", "XGL713_177_sorted.bam", "XGL713_179_sorted.bam", "XGL713_180_sorted.bam", "XGL713_181_sorted.bam", "XGL713_232_sorted.bam", "XGUAE_124_sorted.bam", "XGUAE_36_sorted.bam", "XGUAE_42_sorted.bam", "XGUAE_43_sorted.bam", "XGUAE_44_sorted.bam", "XGUAE_59_sorted.bam", "XGUAE_65_sorted.bam", "XGUAE_70_sorted.bam", "XGUAE_71_sorted.bam", "XGUAE_72_sorted.bam", "XGUAE_92_sorted.bam", "XGUAE_93_sorted.bam", "XGUAE_97_sorted.bam", "XLJONK_14_sorted.bam", "XL_CPT1_sorted.bam", "XL_CPT2_sorted.bam", "XL_CPT3_sorted.bam", "XL_CPT4_sorted.bam", "XM_1_sorted.bam")

d<-ggplot(data=tab, aes(x=EV1,y=-EV2, label = samp.fieldid, color = samp.color)) +
    # label axis 
    labs(x=expression("Eigenvector 1"), y=expression("-Eigenvector 2")) +
    # legend details
    scale_colour_manual(name="Species", values = c("gray"="gray", "green"="green","purple"="purple","red"="red", "blue"="blue","black" = "black"),breaks=c("gray", "green","purple","red","blue","black"),labels=c("X. victorianus", "X. laevis Beaufort West", "X. laevis Niewoudtville","X. gilli","X. laevis Cape Region","X. muelleri"))+
    # add points and fieldID labels
    geom_text_repel(aes(EV1,-EV2, label=(samp.fieldid))) + geom_point(size=4) + 
    # change to cleaner theme
    theme_classic(base_size = 16) +
    # make it clean
    theme_bw()+ theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
    # italicize species names
    theme(legend.text = element_text(face="italic"))+ 
    # move the legend
    theme(legend.position = c(.14, .28)) +
    # add a title
    ggtitle("Principal Components Analsis") + 
    # remove boxes around legend symbols
    theme(legend.key = element_blank())
    #+
    #annotate(geom = "text", x = .2, y = -.4, label = "Southwest Sulawesi", color = "black", angle = -60)+
    #annotate(geom = "text", x = .16, y = .22, label = "North Sulawesi", color = "black", angle = 55)+
    #annotate(geom = "text", x = .24, y = .03, label = "Southeast\n    Sulawesi", color = "black", angle = 20)+
    #annotate(geom = "text", x = -.25, y = -.08, label = "Borneo", color = "black", angle = 10)+
    #annotate(geom = "text", x = -.3, y = .04, label = "Sumatra/Mentawai", color = "black", angle = 10)+
    #annotate(geom = "text", x = -.57, y = 0, label = "Peninsular Malaysia", color = "black", angle = 0)
d
dev.off()

```

# Exclude XM to see what PCA looks like with only XL, XV, and XG

```
~/jre1.8.0_111/bin/java -Xmx2g -jar /home/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T SelectVariants -R ../../XL_v9.1/Xla.v91.repeatMasked.fa --exclude_sample_file XM_sample -o XLXVXG_merged_sorted.bam.vcf.gz --variant XLXVXGXM_merged_sorted.bam.vcf.gz 
```

