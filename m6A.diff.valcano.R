###EXOMEPEAK -> m6A diff figure

python dict_array.py diff_peak.xls  > tmp

awk 'NR==FNR{s[$1]=$0;next}NR>FNR{print s[$1]"\t"$2"\t"$3}' ~/project/project/00.DATABASE/hg38/gene_id_name.3.txt tmp  > tmp.txt

library(ggrepel)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
theme_set(theme_cowplot())

data = read.table("m6A.diff.txt",head=T)
d1d = subset(data,data$fc > 1 & abs(data$fdr) > log(0.05,10), select=c(name,fdr,fc))
d1d$Regulation = 'Upregulated'
d1p = subset(data,data$fc < -1 & abs(data$fdr) > log(0.05,10) , select=c(name,fdr,fc))
d1p$Regulation = 'Downregulated'
d1n = subset(data,abs(data$fdr) < -log(0.05,10) | abs(data$fc) < 1 , select=c(name,fdr,fc))
d1n$Regulation = 'Insignificant'
d1 = rbind(d1d,d1p,d1n)
d1$fdr = abs(d1$fdr)
pdf('m6A.diff.pdf', w=9.5, h= 8, pointsize= 10)
ggplot(d1,aes(fc,fdr,col=Regulation)) + geom_point(size=1,shape=1,alpha=.5) + 
#  facet_grid(~cluster) +
  scale_color_manual(values=c("#377EB8",'grey','#E41A1C')) + xlim(-15,15) +ylim(0,400)+
#  geom_hline(yintercept = -log10(0.05),linetype="dashed") +
  geom_vline(xintercept = 0) + theme_cowplot() + xlab("avg log FC (patient/normal)") +
  xlab("log 2 (fold change)")+
  ylab("-log 10 (adjust P)")
dev.off()

