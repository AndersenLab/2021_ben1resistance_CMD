library(tidyverse)
library(cowplot)
library(ggpubr)
library(glue)

##FIGURE 1 PLOT 
colsnewben1 <- c("N2" = "orange","ECA882"="grey","919"="red","920"="red","1325"="yellow","1326"="yellow","1327"="purple","1328"="purple","1082" = "blue","1081" = "blue", "1076" = "green","1075" = "green", "Deletion" = "grey", "F200Y" = "Red", "E198V"= "Purple","E198L"="yellow","E198A"="blue","F167Y"="green","ECA2795"="cadetblue3","ECA2799"="cadetblue3","ECA2798"="tan","ECA2805"="tan","ECA2804"="pink","ECA2801"="pink","ECA2816"="darkred","ECA2817"="darkred")

##FBZ
dirs <- c("~/Desktop/20210408_newben1/20210409_ABZD2d/", "~/Desktop/20210408_newben1/20210409_ABZE2e/","~/Desktop/20210408_newben1/20210409_ABZF2f/")

full <- rio::import("~/Desktop/2021_newben1paper/data/S1.tsv")

statsfbzd <- full %>%
  dplyr::filter(condition == "fbz30")%>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group2 == 'N2')


medianfbz <- full %>%
  dplyr::filter(trait == "median.EXT")%>%
  #dplyr::filter(condition == "fbz30")%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "ECA882","ECA2795","ECA2799","ECA2798","ECA2805","ECA2801","ECA2804","ECA2816","ECA2817")))%>%
  ggplot+
  aes(x=fancy_strain, y=phenotype)+
  geom_jitter(width=0.1, size = 1)+
  geom_boxplot(aes(fill = fancy_strain, alpha = 0.4),outlier.shape = NA)+
  ylab("Regressed median OD")+
  ggtitle("Fenbendazole")+
  scale_x_discrete(labels=c("N2" = "WT", "ECA882" = "Deletion","ECA2795" = "E198I","ECA2799"="E198I","ECA2798"="E198K","ECA2805"="E198K","ECA2801"="E198T","ECA2804"="E198T","ECA2816"="E198*","ECA2817"="E198*"))+
  stat_pvalue_manual(statsfbzd, label = "p.adj.signif",y.position = c(200),xmax="group1", remove.bracket = TRUE)+
  scale_fill_manual(name = "fancy_strain", labels = c("N2" = "Suceptible","ECA882"="Deletion","ECA2795"="E198I","ECA2798"="E198K", "ECA2799" = "E198I", "ECA2801" = "E198T","ECA2804" = "E198T", "ECA2805" = "E198K","ECA2816" = "E198*", "ECA2817" = "E198*"), values = colsnewben1)+
  cowplot::theme_cowplot(12)+
  theme(legend.position = "None",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())

#F1B
fulla <- rio::import("~/Desktop/2021_newben1paper/data/S2.tsv")

statsabzd <- fulla %>%
  dplyr::filter(trait == "median.EXT")%>%
  dplyr::filter(condition == "abz30")%>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group2 == "N2")

medianabz <- fulla %>%
  dplyr::filter(trait == "median.EXT")%>%
  dplyr::filter(condition == "abz30")%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "ECA882","ECA2795","ECA2799","ECA2798","ECA2805","ECA2801","ECA2804","ECA2816","ECA2817")))%>%
  ggplot+
  aes(x=fancy_strain, y=phenotype)+
  geom_jitter(width=0.1,size=1)+
  geom_boxplot(aes(fill = fancy_strain, alpha = 0.4),outlier.shape = NA)+
  xlab("Allele")+
  ylab("Regressed median OD")+
  ggtitle("Albendazole")+
  stat_pvalue_manual(statsabzd, label = "p.adj.signif",y.position = c(200),xmax="group1", remove.bracket = TRUE)+
  scale_x_discrete(labels=c("N2" = "WT", "ECA882" = "Deletion","ECA2795" = "E198I","ECA2799"="E198I","ECA2798"="E198K","ECA2805"="E198K","ECA2801"="E198T","ECA2804"="E198T","ECA2816"="E198*","ECA2817"="E198*"))+
  scale_fill_manual(name = "fancy_strain", labels = c("N2" = "Suceptible","ECA882"="Deletion","ECA2795"="E198I","ECA2798"="E198K", "ECA2799" = "E198I", "ECA2801" = "E198T","ECA2804" = "E198T", "ECA2805" = "E198K","ECA2816" = "E198*", "ECA2817" = "E198*"), values = colsnewben1)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","ECA882"="Deletion","ECA2795"="E198I","ECA2798"="E198K", "ECA2799" = "E198I", "ECA2801" = "E198T","ECA2804" = "E198T", "ECA2805" = "E198K","ECA1816" = "E198*", "ECA2817" = "E198*"), values = colsnewben1)+
  cowplot::theme_cowplot(12)+
  theme(legend.position = "None",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())

plota <- cowplot::plot_grid(medianabz,medianfbz, labels = c("A","B"))

ggsave("~/Desktop/2021_newben1paper/plots/figure1.png", plot = plota, width = 7.5, height = 3, units = "in")
##FIGURE 3 PLOT
non_contamgfp <- rio::import("~/Desktop/2021_newben1paper/data/S3.tsv")

density_gfp <- density(non_contamgfp$norm.green)

cutoff <- optimize(approxfun(density_gfp$x,density_gfp$y),interval=c(0,2))$minimum

gfp_noncontam <- non_contamgfp %>%
  dplyr::mutate(strain1 = case_when((strain == "N2") ~ "N2",
                                    (strain == "N2") ~ "N2",
                                    (strain == "ECA882") ~ "deletion",
                                    (strain == "ECA882") ~ "deletion",
                                    (strain == "ECA919") ~ "F200Y",
                                    (strain == "ECA919") ~ "F200Y",
                                    (strain == "ECA1076") ~ "F167Y",
                                    (strain == "ECA1076") ~ "F167Y",
                                    (strain == "ECA1082") ~ "E198A",
                                    (strain == "ECA1082") ~ "E198A",
                                    (strain == "ECA1325") ~ "E198L",
                                    (strain == "ECA1325") ~ "E198L",
                                    (strain == "ECA1327") ~ "E198V",
                                    (strain == "ECA1327") ~ "E198V"))%>%
  dplyr::mutate(strain = case_when((strain == "N2" & norm.green > cutoff) ~ "N2(Het)",
                                   (strain == "N2" & norm.green < cutoff) ~ "N2(Hom)",
                                   (strain == "ECA882" & norm.green > cutoff) ~ "del(Het)",
                                   (strain == "ECA882" & norm.green < cutoff) ~ "del(Hom)",
                                   (strain == "ECA919" & norm.green > cutoff) ~ "F200Y(Het)",
                                   (strain == "ECA919" & norm.green < cutoff) ~ "F200Y(Hom)",
                                   (strain == "ECA1076" & norm.green > cutoff) ~ "F167Y(Het)",
                                   (strain == "ECA1076" & norm.green < cutoff) ~ "F167Y(Hom)",
                                   (strain == "ECA1082" & norm.green > cutoff) ~ "E198A(Het)",
                                   (strain == "ECA1082" & norm.green < cutoff) ~ "E198A(Hom)",
                                   (strain == "ECA1325" & norm.green > cutoff) ~ "E198L(Het)",
                                   (strain == "ECA1325" & norm.green < cutoff) ~ "E198L(Hom)",
                                   (strain == "ECA1327" & norm.green > cutoff) ~ "E198V(Het)",
                                   (strain == "ECA1327" & norm.green < cutoff) ~ "E198V(Hom)"))%>%
  dplyr::filter(!(is.na(strain1)))

summedplategfp <- easysorter::sumplate(gfp_noncontam, v3_assay = TRUE, quantiles = TRUE)

bioprunedgfp <- easysorter::bioprune(summedplategfp)

out_prunedgfp <- easysorter::prune_outliers(bioprunedgfp, iqr = FALSE)

fullgfp <- easysorter::regress(out_prunedgfp, assay = FALSE)

fullstatstablegfp <- fullgfp %>%
  dplyr::filter(trait == "median.EXT")%>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd()

write_csv(fullstatstablegfp, file = "~/Desktop/2021_newben1paper/supplemental/supplemental_table7.tsv")

cols <- c("N2(Het)"="orange","N2(Hom)"="orange","del(Het)"="grey","del(Hom)"="grey","F200Y(Het)"="red","F200Y(Hom)"="red","F167Y(Het)"="green","F167Y(Hom)"="green","E198A(Het)"="blue","E198A(Hom)"="blue","E198V(Het)"="purple","E198V(Hom)"="purple","E198L(Het)"="yellow","E198L(Hom)"="yellow")

plotstats <- fullgfp %>%
  dplyr::filter(trait == "median.EXT")%>%
  aov(phenotype ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter((group1 == "N2(Het)" & group2 == "N2(Hom)") | (group1 == "del(Het)" & group2 == "del(Hom)") | (group1 == "F200Y(Het)" & group2 == "F200Y(Hom)") | (group1 == "E198A(Het)" & group2 == "E198A(Hom)") | (group1 == "E198L(Het)" & group2 == "E198L(Hom)")|(group1 == "E198V(Het)" & group2 == "E198V(Hom)")|(group1 == "F167Y(Het)" & group2 == "F167Y(Hom)"))

boxplot <- fullgfp%>%
  dplyr::filter(trait == "median.EXT")%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2(Hom)","N2(Het)","del(Hom)","del(Het)","F200Y(Hom)","F200Y(Het)","F167Y(Hom)","F167Y(Het)","E198A(Hom)","E198A(Het)","E198V(Hom)","E198V(Het)","E198L(Hom)","E198L(Het)")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  scale_x_discrete(breaks = c("N2(Hom)","N2(Het)","del(Hom)","del(Het)","F200Y(Hom)","F200Y(Het)","F167Y(Hom)","F167Y(Het)","E198A(Hom)","E198A(Het)","E198V(Hom)","E198V(Het)","E198L(Hom)","E198L(Het)"),labels = c("+\n+","+\n+","Δben-1\nΔben-1","Δben-1\n+","F200Y\nF200Y","F200Y\n+","F167Y\nF167Y","F167Y\n+","E198A\nE198A","E198A\n+","E198V\nE198V","E198V\n+","E198L\nE198L","E198L\n+"))+
  geom_jitter(width = 0.1,size =1)+
  geom_boxplot(aes(fill=fancy_strain, alpha = 0.4),outlier.shape = NA)+
  stat_pvalue_manual(plotstats, label = "p.adj.signif",xmin="group1",xmax="group2", y.position = c(550),remove.bracket = TRUE,size = 4)+
  scale_fill_manual(name = "fancy_strain",values = cols)+
  ylab("Regressed median OD")+
  cowplot::theme_cowplot(12)+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position = "none")

GFPsplit <- non_contamgfp %>%
  dplyr::mutate(GFP = ifelse(non_contamgfp$norm.green > cutoff,"GFP","nonGFP"))%>%
  ggplot()+
  aes(x=norm.green)+
  geom_density(aes(fill=GFP))+
  scale_fill_manual(values = c("GFP"="chartreuse4", "nonGFP"="grey"))+
  geom_vline(xintercept = cutoff)+
  xlab("Normalized GFP expression")+
  ylab("")+
  ylim(0,8)+
  geom_text(x=0.25,y=8,label="GFP(-)")+
  geom_text(aes(fontface = "italic",x=0.25,y=7),label="ben-1")+
  geom_text(aes(fontface ="italic",x=0.25,y=6.3),label="ben-1")+
  geom_segment(aes(x = -0.1, y = 6.6, xend = 0.5, yend = 6.6))+
  geom_text(aes(x=2.25,y=3.5),label="GFP(+)")+
  geom_text(aes(x=2.25,y=2.5, fontface = "italic"),label="ben-1")+
  geom_text(x=2.25,y=1.8,label="+")+
  geom_segment(aes(x = 1.9, y = 2.1, xend = 2.5, yend = 2.1))+
  cowplot::theme_cowplot(12)+
  theme(legend.position = "none")

scheme <- ggdraw()+draw_image("~/Desktop/2021_newben1paper/plots/Crosssetup.png")

fig3 <- cowplot::plot_grid(cowplot::plot_grid(scheme,GFPsplit,ncol = 2,labels = c("A","B")), boxplot,nrow = 2,labels = c("A","C"))

ggsave("~/Desktop/2021_newben1paper/plots/figure3.png", plot = fig3,width = 8.5,height = 6,units = "in")

##Supplemental Figure 1
fullall <- rio::import("~/Desktop/2021_newben1paper/data/S4.tsv")

fullallstatsabz <- fullall %>%
  dplyr::filter(condition == "abz30")%>%
  dplyr::filter(trait == "median.EXT")%>%
  aov(phenotype ~ strain, data=.)%>%
  rstatix::tukey_hsd()%>%
  dplyr::mutate(condition = "Albendazole")

fullallstatsfbz <- fullall %>%
  dplyr::filter(condition == "fbz30")%>%
  dplyr::filter(trait == "median.EXT")%>%
  aov(phenotype ~ strain, data=.)%>%
  rstatix::tukey_hsd()%>%
  dplyr::mutate(condition = "Fenbendazole")

allstats <- rbind(fullallstatsabz,fullallstatsfbz)
write_csv(allstats, file = "~/Desktop/2021_newben1paper/supplemental/supplemental_table4.tsv")
fullabzN2 <-fullallstatsabz %>%
  dplyr::filter(group2 == "N2")

fullfbzN2 <- fullallstatsfbz %>%
  dplyr::filter(group2 == "N2")

allallelesabzconditions <- fullall%>%
  dplyr::filter(trait == "median.EXT")%>%
  dplyr::filter(condition == "abz30")%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "ECA882","ECA2795","ECA2799","ECA2798","ECA2805","ECA2801","ECA2804","ECA2816","ECA2817")))%>%
  ggplot+
  aes(x=fancy_strain, y=phenotype)+
  ggtitle("Albendazole")+
  geom_jitter(width=0.1, size = 1)+
  geom_boxplot(aes(fill = fancy_strain, alpha = 0.4),outlier.shape = NA)+
  ylab("Regressed median OD")+
  stat_pvalue_manual(fullabzN2,size = 4, label = "p.adj.signif",y.position = c(200),xmax="group1", remove.bracket = TRUE)+
  scale_x_discrete(labels=c("N2" = "WT", "ECA882" = "ean64","ECA2795" = "ean217","ECA2799"="ean219","ECA2798"="ean218","ECA2805"="ean213","ECA2801"="ean210","ECA2804"="ean209","ECA2816"="ean214","ECA2817"="ean215"))+
  scale_fill_manual(name = "fancy_strain", labels = c("N2" = "Suceptible","ECA882"="Deletion","ECA2795"="E198I","ECA2798"="E198K", "ECA2799" = "E198I", "ECA2801" = "E198T","ECA2804" = "E198T", "ECA2805" = "E198K","ECA2816" = "E198*", "ECA2817" = "E198*"), values = colsnewben1)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","ECA882"="Deletion","ECA2795"="E198I","ECA2798"="E198K", "ECA2799" = "E198I", "ECA2801" = "E198T","ECA2804" = "E198T", "ECA2805" = "E198K","ECA1816" = "E198*", "ECA2817" = "E198*"), values = colsnewben1)+
  cowplot::theme_cowplot(12)+
  theme(legend.position = "None",
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())

allallelesfbzconditions <- fullall%>%
  dplyr::filter(trait == "median.EXT")%>%
  dplyr::filter(condition == "fbz30")%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "ECA882","ECA2795","ECA2799","ECA2798","ECA2805","ECA2801","ECA2804","ECA2816","ECA2817")))%>%
  ggplot+
  aes(x=fancy_strain, y=phenotype)+
  geom_jitter(width=0.1, size = 1)+
  geom_boxplot(aes(fill = fancy_strain, alpha = 0.4),outlier.shape = NA)+
  ylab("Regressed median OD")+
  ggtitle("Fenbendazole")+
  stat_pvalue_manual(fullfbzN2,size=4, label = "p.adj.signif",y.position = c(200),xmax="group1", remove.bracket = TRUE)+
  scale_x_discrete(labels=c("N2" = "WT", "ECA882" = "ean64","ECA2795" = "ean217","ECA2799"="ean219","ECA2798"="ean218","ECA2805"="ean213","ECA2801"="ean210","ECA2804"="ean209","ECA2816"="ean214","ECA2817"="ean215"))+
  scale_fill_manual(name = "fancy_strain", labels = c("N2" = "Suceptible","ECA882"="Deletion","ECA2795"="E198I","ECA2798"="E198K", "ECA2799" = "E198I", "ECA2801" = "E198T","ECA2804" = "E198T", "ECA2805" = "E198K","ECA2816" = "E198*", "ECA2817" = "E198*"), values = colsnewben1)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","ECA882"="Deletion","ECA2795"="E198I","ECA2798"="E198K", "ECA2799" = "E198I", "ECA2801" = "E198T","ECA2804" = "E198T", "ECA2805" = "E198K","ECA1816" = "E198*", "ECA2817" = "E198*"), values = colsnewben1)+
  cowplot::theme_cowplot(12)+
  theme(legend.position = "None",
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())


allalleles <- cowplot::plot_grid(allallelesabzconditions, allallelesfbzconditions, ncol=2, labels = c("A","B"))
ggsave("~/Desktop/2021_newben1paper/plots/supplementalfigure1.png", plot = allalleles, width = 7.5,height = 6,units = "in")

##Supplemental figure 2

outprunedDMSO <- rio::import("~/Desktop/2021_newben1paper/data/S5.tsv")

fullallstatsDMSO <- outprunedDMSO %>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(trait == "median.EXT")%>%
  aov(phenotype ~ strain, data=.)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group2=="N2")


DMSOplot <- outprunedDMSO%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(trait == "median.EXT")%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "ECA882","ECA2795","ECA2799","ECA2798","ECA2805","ECA2801","ECA2804","ECA2816","ECA2817")))%>%
  ggplot()+
  aes(x=fancy_strain, y=phenotype)+
  geom_jitter(width=0.1)+
  geom_boxplot(aes(fill = fancy_strain, alpha=0.4),outlier.shape = NA)+
  ylab("Regressed median OD")+
  stat_pvalue_manual(fullallstatsDMSO, label = "p.adj.signif",y.position = c(140),xmax="group1", remove.bracket = TRUE)+
  scale_x_discrete(labels=c("N2" = "WT", "ECA882" = "Deletion","ECA2795" = "E198I","ECA2799"="E198I","ECA2798"="E198K","ECA2805"="E198K","ECA2801"="E198T","ECA2804"="E198T","ECA2816"="E198*","ECA2817"="E198*"))+
  scale_fill_manual(name = "fancy_strain", labels = c("N2" = "Suceptible","ECA882"="Deletion","ECA2795"="E198I","ECA2798"="E198K", "ECA2799" = "E198I", "ECA2801" = "E198T","ECA2804" = "E198T", "ECA2805" = "E198K","ECA2816" = "E198*", "ECA2817" = "E198*"), values = colsnewben1)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","ECA882"="Deletion","ECA2795"="E198I","ECA2798"="E198K", "ECA2799" = "E198I", "ECA2801" = "E198T","ECA2804" = "E198T", "ECA2805" = "E198K","ECA1816" = "E198*", "ECA2817" = "E198*"), values = colsnewben1)+
  cowplot::theme_cowplot(12)+
  theme(legend.position = "None",
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank())
ggsave("~/Desktop/2021_newben1paper/plots/supplementalfigure2.png", plot = DMSOplot, width = 7.5, height = 3, units = "in")

##Supplemental figure 3


fullallstatsDMSOgfp <- out_prunedgfp %>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(trait == "median.EXT")%>%
  aov(phenotype ~ strain, data=.)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group2 %in% c("N2(Hom)","N2(Het)"))%>%
  dplyr::filter(group1 != "N2(Het)")


DMSOGFP <- out_prunedgfp%>%
  dplyr::filter(trait == "median.EXT")%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2(Hom)","N2(Het)","del(Hom)","del(Het)","F200Y(Hom)","F200Y(Het)","F167Y(Hom)","F167Y(Het)","E198A(Hom)","E198A(Het)","E198V(Hom)","E198V(Het)","E198L(Hom)","E198L(Het)")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  scale_x_discrete(breaks = c("N2(Hom)","N2(Het)","del(Hom)","del(Het)","F200Y(Hom)","F200Y(Het)","F167Y(Hom)","F167Y(Het)","E198A(Hom)","E198A(Het)","E198V(Hom)","E198V(Het)","E198L(Hom)","E198L(Het)"),labels = c("+\n+","+\n+","Δben-1\nΔben-1","Δben-1\n+","F200Y\nF200Y","F200Y\n+","F167Y\nF167Y","F167Y\n+","E198A\nE198A","E198A\n+","E198V\nE198V","E198V\n+","E198L\nE198L","E198L\n+"))+
  geom_jitter(width = 0.1,size =1)+
  geom_boxplot(aes(fill=fancy_strain, alpha = 0.4),outlier.shape = NA)+
  stat_pvalue_manual(fullallstatsDMSOgfp, label = "p.adj.signif",y.position = c(1600),xmax="group1", remove.bracket = TRUE)+
  scale_fill_manual(name = "fancy_strain",values = cols)+
  ylab("Normalized median OD")+
  cowplot::theme_cowplot(12)+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position = "none")

ggsave("~/Desktop/2021_newben1paper/plots/supplementalfigure3.png", plot=DMSOGFP, width = 8.5, height = 4, units = "in")


stats <- NULL
for (i in unique(gfp_noncontam$strain)) {
  x <- gfp_noncontam %>%
    dplyr::filter(strain == i)%>%
    aov(TOF ~ condition, data = .)
  xfull <- as.data.frame(TukeyHSD(x)[[1]])
  xfull <- rownames_to_column(xfull, var = "comp")%>%
    dplyr::mutate(strain = i)
  stats <- rbind(stats, xfull)
  rm(x,xfull)
}

supplementalfigure4 <- gfp_noncontam%>%
  dplyr::filter(!(is.na(condition)))%>%
  dplyr::filter(!(is.na(strain)))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2(Hom)","N2(Het)","del(Hom)","del(Het)","F200Y(Hom)","F200Y(Het)","F167Y(Hom)","F167Y(Het)","E198A(Hom)","E198A(Het)","E198V(Hom)","E198V(Het)","E198L(Hom)","E198L(Het)")))%>%
  dplyr::mutate(statsig = case_when(strain == "N2(Hom)" ~ "NS",
                                    strain == "N2(Het)" ~ "NS",
                                    strain == "del(Hom)" ~ "NS",
                                    strain == "del(Het)" ~ "****",
                                    strain == "F200Y(Hom)" ~"NS",
                                    strain == "F200Y(Het)" ~ "****",
                                    strain == "F167Y(Hom)" ~ "NS",
                                    strain == "F167Y(Het)" ~ "****",
                                    strain == "E198A(Hom)" ~ "NS",
                                    strain == "E198A(Het)" ~ "****",
                                    strain == "E198L(Hom)" ~ "NS",
                                    strain == "E198L(Het)" ~"****",
                                    strain == "E198V(Hom)" ~  "NS",
                                    strain == "E198V(Het)" ~ "****"))%>%
  dplyr::mutate(Condition = case_when(condition == "DMSO" ~ "Control",
                                      condition == "abz30" ~ "Albendazole"))%>%
  ggplot()+
  aes(x=fancy_strain,y = TOF)+
  geom_boxplot(aes(fill = Condition), outlier.shape = NA)+
  facet_grid( ~ fancy_strain, scales = "free")+
  cowplot::theme_cowplot(12)+
  ylim(100,750)+
  scale_x_discrete(breaks = c("N2(Hom)","N2(Het)","del(Hom)","del(Het)","F200Y(Hom)","F200Y(Het)","F167Y(Hom)","F167Y(Het)","E198A(Hom)","E198A(Het)","E198V(Hom)","E198V(Het)","E198L(Hom)","E198L(Het)"),labels = c("+\n+","+\n+","Δben-1\nΔben-1","Δben-1\n+","F200Y\nF200Y","F200Y\n+","F167Y\nF167Y","F167Y\n+","E198A\nE198A","E198A\n+","E198V\nE198V","E198V\n+","E198L\nE198L","E198L\n+"))+
  ylab("Worm legnth")+
  xlab("Strain")+
  geom_text(aes(x=fancy_strain, y=725, label = statsig))+
  theme(legend.position = "top",
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank())
ggsave("~/Desktop/2021_newben1paper/plots/supplementalfigure4.png", plot = supplementalfigure4, width = 7.5, height = 4, units = "in")

