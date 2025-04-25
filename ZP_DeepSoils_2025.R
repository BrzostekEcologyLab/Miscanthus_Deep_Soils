#### Packages ####

library(tidyverse)
library(ggplot2)
library(multcompView)
library(cowplot)
library(patchwork)
library(stringr)
library(grid) #table for LME
library(gridExtra) #table grob for LME
library(performance) #checking for collinearity
library(nlme)
library(lme4)
library(marginaleffects)
library(car)
library(MuMIn)
library(broom)

#### Functions and Themes ####

SE.1 <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  (sqrt(var(x)/length(x)))
}

Mean.SE <- list(mean, SE.1)
names(Mean.SE) <- c("Mean", "SE")

theme_zp <- function() {
  theme_classic() +
    theme(legend.position = "none",
          axis.title.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12))
}

#### Elements ####

zp <- read.csv("Pagliaro_DeepSoils_Data.csv")

zp <- zp %>% 
  mutate(across(c(pit, depth), as.factor))
str(zp)

zp$depth <- ordered(zp$depth, levels = c("0 to 15", "15 to 30", "30 to 50", 
                                         "50 to 70", "70 to 100"))
  #making depth an ordered factor for downstream analysis

zp$vol_cm3 <- as.numeric(zp$vol_cm3)

zpPalette <- c("#808080") #monotone
zpPalette2 <- c("#808080","#561F37") #duotone


zp.sum <- zp %>% #for ggplot2
  group_by(depth) %>%
  reframe(
    across(AP:MAOC_change_mmol, Mean.SE, .unpack = TRUE, na.rm = TRUE)
  )

#### Statistics ####

vars <- colnames(zp[,9:36])

zp.fp <- as.data.frame(matrix(,28,2))
rownames(zp.fp) <- vars
colnames(zp.fp) <- c("F-stat", "p-val")

options(na.action = "na.omit") #for the anovas (changed below for LME model)

  for(i in 1:28) {
    zp.fp[i,] <- summary(aov(
      get(vars[i]) ~ depth, zp), na.rm = TRUE)[[1]][1,4:5] #F and p table
  }

zp.res <- vector("list", length = 28) #list for all anovas

for(i in 1:28) {
  zp.res[[i]][[1]] <- list(aov(get(vars[i]) ~ depth, zp)) #all anovas in a list
}

zp.lt <- rownames_to_column(as.data.frame(
  multcompLetters4(zp.res[[1]][[1]][[1]], 
                   TukeyHSD(zp.res[[1]][[1]][[1]]))$depth$Letters))


for(i in 2:28) {
  zp.lt <- left_join(zp.lt, (rownames_to_column(as.data.frame(
    multcompLetters4(zp.res[[i]][[1]][[1]], 
                     TukeyHSD(zp.res[[i]][[1]][[1]]))$depth$Letters))), 
    by = "rowname")
}

zp.lt <- column_to_rownames(zp.lt, var = "rowname")

colnames(zp.lt) <- paste0(vars, "_Let")
#rownames(zp.lt) <- zp.sum$depth
zp.lt$depth <- rownames(zp.lt)

zp.gg <- full_join(zp.sum, zp.lt, by = c("depth"))

zp.gg <- zp.gg %>% 
  rename_with(~str_replace(., '_Mean', ''))

#### Letters ####
#Not really needed since letters are autopopulated & added to zp.gg! BUT:

a <- aov(gC_cm3soil ~ depth, zp) #Total C (Fig 1 panel A)
a.tk <- TukeyHSD(a)
a.let <- multcompLetters4(a, a.tk)
a.let$depth #0 to 100: a, a, b, bc, c
rm(a, a.tk, a.let)

a <- aov(gMAOC_cm3soil ~ depth, zp) #MAOC (Fig 1 panel B)
a.tk <- TukeyHSD(a)
a.let <- multcompLetters4(a, a.tk)
a.let$depth #0 to 100: ab, a, abc, bc, c
rm(a, a.tk, a.let)

a <- aov(gPOC_cm3soil ~ depth, zp) #POC (Fig 1 panel C)
a.tk <- TukeyHSD(a)
a.let <- multcompLetters4(a, a.tk)
a.let$depth #0 to 100: a, ab, ab, b, b
rm(a, a.tk, a.let)

a <- aov(MAOC_POC ~ depth, zp) #MAOC:POC (Fig 1 panel D)
a.tk <- TukeyHSD(a)
a.let <- multcompLetters4(a, a.tk)
a.let$depth #0 to 100: b, ab, ab, a, a
rm(a, a.tk, a.let)

a <- aov(froot_g_cm3soil ~ depth, zp) #Fine root biomass (Table 1 row 1)
a.tk <- TukeyHSD(a)
a.let <- multcompLetters4(a, a.tk)
a.let$depth #0 to 100: a, b, b, b, b
rm(a, a.tk, a.let)

#### Table 1 ####

table.1 <- select(zp.gg, depth, 
                  froot_g_cm3soil, froot_g_cm3soil_SE, froot_g_cm3soil_Let,
                  R_mmol_nosub, R_mmol_nosub_SE, R_mmol_nosub_Let,
                  Min_ug_gsoil_day, Min_ug_gsoil_day_SE, Min_ug_gsoil_day_Let,
                  Nit_ug_gsoil_day, Nit_ug_gsoil_day_SE, Nit_ug_gsoil_day_Let,
                  NAG, NAG_SE, NAG_Let, 
                  BG, BG_SE, BG_Let, 
                  AP, AP_SE, AP_Let) #df with means, SE, and letter reports

table.1.vars <- c("froot_g_cm3soil","R_mmol_nosub","Min_ug_gsoil_day",
                  "Nit_ug_gsoil_day","NAG","BG","AP")

table.1.fp <- zp.fp %>%
  rownames_to_column() %>%
  filter(rowname %in% table.1.vars) #df with F-stat and p-vals

write.csv(table.1, file = "Table_1_MeanSE.csv", row.names = F) 

#### Graphing ####

APplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=AP))+ #1
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=AP-AP_SE,ymax=AP+AP_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="AP Activity (umol/g/hr)",x="Soil Depth (cm)") +
  geom_text(aes(label = AP_Let, y = AP+AP_SE), position = position_dodge(0.9), 
            hjust = -0.5, size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=0.35, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("AP"),][2], digits = 3))), 
           color="black")+
  coord_flip()
APplot

ggsave("APplot.png", plot = APplot, dpi = 300, height = 6, width = 5)

BGplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=BG))+ #2
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=BG-BG_SE,ymax=BG+BG_SE), width=0.2, 
                position=position_dodge(0.9)) + 
  theme_zp() +
  scale_fill_manual(values = zpPalette) + 
  labs(y="BG Activity (umol/g/hr)",x="Soil Depth (cm)") +
  geom_text(aes(label = BG_Let, y = BG+BG_SE), position = position_dodge(0.9), 
            hjust = -0.5, size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=0.35, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("BG"),][2], digits = 3))), 
           color="black")+
  coord_flip()
BGplot

ggsave("BGplot.png", plot = BGplot, dpi = 300, height = 6, width = 5)

NAGplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=NAG))+ #3
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=NAG-NAG_SE,ymax=NAG+NAG_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() +
  scale_fill_manual(values = zpPalette) + 
  labs(y="NAG Activity (umol/g/hr)",x="Soil Depth (cm)") +
  geom_text(aes(label = NAG_Let, y = NAG+NAG_SE), position = 
              position_dodge(0.9), hjust = -0.5, size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=0.020, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("NAG"),][2], digits = 3))), 
           color="black")+
  coord_flip()
NAGplot

ggsave("NAGplot.png", plot = NAGplot, dpi = 300, height = 6, width = 5)

Tot_C_gplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=Tot_C_g))+ #4
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=Tot_C_g-Tot_C_g_SE,ymax=Tot_C_g+Tot_C_g_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() +
  scale_fill_manual(values = zpPalette) + 
  labs(y="Carbon (g·soil g-)",x="Soil Depth (cm)") +
  geom_text(aes(label = Tot_C_g_Let, y = Tot_C_g+Tot_C_g_SE), 
            position = position_dodge(0.9), hjust = -0.5, 
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=1500, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("Tot_C_g"),][2], digits = 3))), 
           color="black")+
  coord_flip()
Tot_C_gplot

ggsave("Tot_C_gplot.png", plot = Tot_C_gplot, dpi = 300, height = 6, width = 5)

gC_cm3soilplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=gC_cm3soil))+ #5
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=gC_cm3soil-gC_cm3soil_SE,ymax=gC_cm3soil+gC_cm3soil_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() +
  scale_fill_manual(values = zpPalette) + 
  labs(y="Carbon (g·soil cm3)",x="Soil Depth (cm)") +
  geom_text(aes(label = gC_cm3soil_Let, y = gC_cm3soil+gC_cm3soil_SE), 
            position = position_dodge(0.9), hjust = -0.5, 
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=0.035, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("gC_cm3soil"),][2], digits = 3))), 
           color="black")+
  coord_flip()
gC_cm3soilplot

ggsave("gC_cm3soilplot.png", plot = gC_cm3soilplot, dpi = 300, height = 6, 
       width = 5)

Tot_N_gplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=Tot_N_g))+ #6
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=Tot_N_g-Tot_N_g_SE,ymax=Tot_N_g+Tot_N_g_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() +
  scale_fill_manual(values = zpPalette) + 
  labs(y="Nitrogen (g·soil g-)",x="Soil Depth (cm)") +
  geom_text(aes(label = Tot_N_g_Let, y = Tot_N_g+Tot_N_g_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=125, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("Tot_N_g"),][2], digits = 3))), 
           color="black")+
  coord_flip()
Tot_N_gplot

ggsave("Tot_N_gplot.png", plot = Tot_N_gplot, dpi = 300, height = 6, width = 5)

gN_cm3soilplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=gN_cm3soil))+ #7
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=gN_cm3soil-gN_cm3soil_SE,ymax=gN_cm3soil+gN_cm3soil_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() +
  scale_fill_manual(values = zpPalette) + 
  labs(y="Nitrogen (g·soil cm3)",x="Soil Depth (cm)") +
  geom_text(aes(label = gN_cm3soil_Let, y = gN_cm3soil+gN_cm3soil_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=0.0025, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("gN_cm3soil"),][2], digits = 3))), 
           color="black")+
  coord_flip()
gN_cm3soilplot

ggsave("gN_cm3soilplot.png", plot = gN_cm3soilplot, dpi = 300, height = 6, 
       width = 5)

soil_CNplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=soil_CN))+ #8
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=soil_CN-soil_CN_SE,ymax=soil_CN+soil_CN_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() +
  scale_fill_manual(values = zpPalette) + 
  labs(y="Soil C:N",x="Soil Depth (cm)") +
  geom_text(aes(label = soil_CN_Let, y = soil_CN+soil_CN_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=12, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("soil_CN"),][2], digits = 3))), 
           color="black")+
  coord_flip()
soil_CNplot

ggsave("soil_CNplot.png", plot = soil_CNplot, dpi = 300, height = 6, width = 5)

Nit_ug_gsoil_dayplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=Nit_ug_gsoil_day))+ #9
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=Nit_ug_gsoil_day-Nit_ug_gsoil_day_SE,
                    ymax=Nit_ug_gsoil_day+Nit_ug_gsoil_day_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Nitrification (ug/soil g/day)",x="Soil Depth (cm)") +
  geom_text(aes(label = Nit_ug_gsoil_day_Let, 
                y = Nit_ug_gsoil_day+Nit_ug_gsoil_day_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=0.75, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("Nit_ug_gsoil_day"),][2], digits = 3))), 
           color="black")+
  coord_flip()
Nit_ug_gsoil_dayplot

ggsave("Nit_ug_gsoil_dayplot.png", plot = Nit_ug_gsoil_dayplot, dpi = 300, 
       height = 6, width = 5)

Min_ug_gsoil_dayplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=Min_ug_gsoil_day))+ #10
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=Min_ug_gsoil_day-Min_ug_gsoil_day_SE,
                    ymax=Min_ug_gsoil_day+Min_ug_gsoil_day_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Mineralization (ug/soil g/day)",x="Soil Depth (cm)") +
  geom_text(aes(label = Min_ug_gsoil_day_Let, 
                y = Min_ug_gsoil_day+Min_ug_gsoil_day_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=0.75, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("Min_ug_gsoil_day"),][2], 
                                       digits = 3))), 
           color="black")+
  coord_flip()
Min_ug_gsoil_dayplot

ggsave("Min_ug_gsoil_dayplot.png", plot = Min_ug_gsoil_dayplot, dpi = 300, 
       height = 6, width = 5)

froot_g_cm3soilplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=froot_g_cm3soil))+ #11
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=froot_g_cm3soil-froot_g_cm3soil_SE,
                    ymax=froot_g_cm3soil+froot_g_cm3soil_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Fine roots (g·soil cm3)",x="Soil Depth (cm)") +
  geom_text(aes(label = froot_g_cm3soil_Let, 
                y = froot_g_cm3soil+froot_g_cm3soil_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=0.0025, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("froot_g_cm3soil"),][2], 
                                       digits = 3))), 
           color="black")+
  coord_flip()
froot_g_cm3soilplot

ggsave("froot_g_cm3soilplot.png", plot = froot_g_cm3soilplot, dpi = 300, 
       height = 6, width = 5)

root_CNplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=root_CN))+ #12
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=root_CN-root_CN_SE,ymax=root_CN+root_CN_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Root C:N",x="Soil Depth (cm)") +
  geom_text(aes(label = root_CN_Let, y = root_CN+root_CN_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=0.7, y=100, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("root_CN"),][2], digits = 3))), 
           color="black")+
  coord_flip()
root_CNplot

ggsave("root_CNplot.png", plot = root_CNplot, dpi = 300, height = 6, width = 5)

Tot_LF_gplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=Tot_LF_g))+ #13
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=Tot_LF_g-Tot_LF_g_SE,ymax=Tot_LF_g+Tot_LF_g_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Total Light Fraction (g)",x="Soil Depth (cm)") +
  geom_text(aes(label = Tot_LF_g_Let, y = Tot_LF_g+Tot_LF_g_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=1600, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("Tot_LF_g"),][2], digits = 3))), 
           color="black")+
  coord_flip()
Tot_LF_gplot

ggsave("Tot_LF_gplot.png", plot = Tot_LF_gplot, dpi = 300, height = 6, width = 5)

Tot_POM_gplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=Tot_POM_g))+ #14
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=Tot_POM_g-Tot_POM_g_SE,ymax=Tot_POM_g+Tot_POM_g_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Total POM (g)",x="Soil Depth (cm)") +
  geom_text(aes(label = Tot_POM_g_Let, y = Tot_POM_g+Tot_POM_g_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=5, y=10000, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("Tot_POM_g"),][2], digits = 3))), 
           color="black")+
  coord_flip()
Tot_POM_gplot

ggsave("Tot_POM_gplot.png", plot = Tot_POM_gplot, dpi = 300, height = 6, width = 5)

Tot_MAOM_gplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=Tot_MAOM_g))+ #15
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=Tot_MAOM_g-Tot_MAOM_g_SE,ymax=Tot_MAOM_g+Tot_MAOM_g_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Total MAOM (g)",x="Soil Depth (cm)") +
  geom_text(aes(label = Tot_MAOM_g_Let, y = Tot_MAOM_g+Tot_MAOM_g_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=5, y=60000, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("Tot_MAOM_g"),][2], digits = 3))), 
           color="black")+
  coord_flip()
Tot_MAOM_gplot

ggsave("Tot_MAOM_gplot.png", plot = Tot_MAOM_gplot, dpi = 300, height = 6, 
       width = 5)

LF_C_gplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=LF_C_g))+
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=LF_C_g-LF_C_g_SE,ymax=LF_C_g+LF_C_g_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Carbon in Light Fraction (g)",x="Soil Depth (cm)") +
  geom_text(aes(label = LF_C_g_Let, y = LF_C_g+LF_C_g_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=380, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("LF_C_g"),][2], digits = 3))), 
           color="black")+
  coord_flip()
LF_C_gplot

ggsave("LF_C_gplot.png", plot = LF_C_gplot, dpi = 300, height = 6, width = 5)

POM_C_gplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=POM_C_g))+
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=POM_C_g-POM_C_g_SE,ymax=POM_C_g+POM_C_g_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Carbon in POM (g)",x="Soil Depth (cm)") +
  geom_text(aes(label = POM_C_g_Let, y = POM_C_g+POM_C_g_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=200, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("POM_C_g"),][2], digits = 3))), 
           color="black")+
  coord_flip()
POM_C_gplot

ggsave("POM_C_gplot.png", plot = POM_C_gplot, dpi = 300, height = 6, width = 5)

MAOM_C_gplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=MAOM_C_g ))+
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=MAOM_C_g-MAOM_C_g_SE,ymax=MAOM_C_g+MAOM_C_g_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Carbon in MAOM (g)",x="Soil Depth (cm)") +
  geom_text(aes(label = MAOM_C_g_Let, y = MAOM_C_g+MAOM_C_g_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=1400, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("MAOM_C_g"),][2], digits = 3))), 
           color="black")+
  coord_flip()
MAOM_C_gplot

ggsave("MAOM_C_gplot.png", plot = MAOM_C_gplot, dpi = 300, height = 6, width = 5)

LHPOM_C_gplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=LHPOM_C_g))+
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=LHPOM_C_g-LHPOM_C_g_SE,ymax=LHPOM_C_g+LHPOM_C_g_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Carbon in Light and Heavy POM (g)",x="Soil Depth (cm)") +
  geom_text(aes(label = LHPOM_C_g_Let, y = LHPOM_C_g+LHPOM_C_g_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=400, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("LHPOM_C_g"),][2], digits = 3))), 
           color="black")+
  coord_flip()
LHPOM_C_gplot

ggsave("LHPOM_C_gplot.png", plot = LHPOM_C_gplot, dpi = 300, height = 6, width = 5)

MAOC_POCplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=MAOC_POC))+
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=MAOC_POC-MAOC_POC_SE,ymax=MAOC_POC+MAOC_POC_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="MAOM-C:POM-C",x="Soil Depth (cm)") +
  geom_text(aes(label = MAOC_POC_Let, y = MAOC_POC+MAOC_POC_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=5, y=35, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("MAOC_POC"),][2], digits = 3))), 
           color="black")+
  coord_flip()
MAOC_POCplot

ggsave("MAOC_POCplot.png", plot = MAOC_POCplot, dpi = 300, height = 6, width = 5)

Frac_C_gplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=Frac_C_g))+
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=Frac_C_g-Frac_C_g_SE,ymax=Frac_C_g+Frac_C_g_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Carbon in All Fractions (g)",x="Soil Depth (cm)") +
  geom_text(aes(label = Frac_C_g_Let, y = Frac_C_g+Frac_C_g_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=1500, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("Frac_C_g"),][2], digits = 3))), 
           color="black")+
  coord_flip()
Frac_C_gplot

ggsave("Frac_C_gplot.png", plot = Frac_C_gplot, dpi = 300, height = 6, width = 5)

POC_cm3_plot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=gPOC_cm3soil ))+
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=gPOC_cm3soil-gPOC_cm3soil_SE,
                    ymax=gPOC_cm3soil+gPOC_cm3soil_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Carbon in POM (g·cm3)",x="Soil Depth (cm)") +
  geom_text(aes(label = gPOC_cm3soil_Let, y = gPOC_cm3soil+gPOC_cm3soil_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=0.01, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("gPOC_cm3soil"),][2], 
                                       digits = 3))), 
           color="black")+
  coord_flip()
POC_cm3_plot

MAOC_cm3_plot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=gMAOC_cm3soil ))+
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=gMAOC_cm3soil-gMAOC_cm3soil_SE,
                    ymax=gMAOC_cm3soil+gMAOC_cm3soil_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Carbon in MAOM (g·cm3)",x="Soil Depth (cm)") +
  geom_text(aes(label = gMAOC_cm3soil_Let, y = gMAOC_cm3soil+gMAOC_cm3soil_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") +
  annotate(geom="text", x=1, y=0.03, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("gMAOC_cm3soil"),][2], 
                                       digits = 3))), 
           color="black")+
  coord_flip()
MAOC_cm3_plot

R_mmol_nosubplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=R_mmol_nosub ))+
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=R_mmol_nosub-R_mmol_nosub_SE,
                    ymax=R_mmol_nosub+R_mmol_nosub_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) + 
  labs(y="Respiration (mmol CO2·dry soil g) - no substrate",x="Soil Depth (cm)") +
  geom_text(aes(label = R_mmol_nosub_Let, y = R_mmol_nosub+R_mmol_nosub_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") + 
  annotate(geom="text", x=1, y=0.013, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("R_mmol_nosub"),][2], 
                                       digits = 3))), 
           color="black")+
  coord_flip()
R_mmol_nosubplot

ggsave("R_mmol_nosubplot.png", plot = R_mmol_nosubplot, dpi = 300, height = 6, width = 5)

R_mmol_subplot <- ggplot(zp.gg,aes(x=fct_rev(depth), y=R_mmol_sub ))+
  geom_bar(position=position_dodge(),stat="identity") + 
  geom_errorbar(aes(ymin=R_mmol_sub-R_mmol_sub_SE,ymax=R_mmol_sub+R_mmol_sub_SE), 
                width=0.2,position=position_dodge(0.9)) + 
  theme_zp() + 
  scale_fill_manual(values = zpPalette) +
  labs(y="Respiration (mmol 13-CO2·dry soil g)",x="Soil Depth (cm)") +
  geom_text(aes(label = R_mmol_sub_Let, y = R_mmol_sub+R_mmol_sub_SE), 
            position = position_dodge(0.9), hjust = -0.5,
            size = 10, size.unit = "pt") +
  annotate(geom="text", x=1, y=0.004, 
           label=as.character(c(signif(zp.fp[row.names(zp.fp) %in% 
                                               c("R_mmol_sub"),][2], 
                                       digits = 3))), 
           color="black")+
  coord_flip()
R_mmol_subplot

ggsave("R_mmol_subplot.png", plot = R_mmol_subplot, dpi = 300, height = 6, 
       width = 5)

#connecting letters report for the following box plots
zp.gg$MAOM_13C_mmol_Let #Fig 3A (for alpha = 0.05)
zp.gg$MAOM_13Cprop_Let #Fig 3B (for alpha = 0.05)
zp.gg$MAOC_change_mmol_Let  #Fig 3C (for alpha = 0.05)
MAOC_change_mmol.aov <- aov(MAOC_change_mmol ~ depth, zp)
MAOC_change_mmol.tk <- TukeyHSD(MAOC_change_mmol.aov) #p-adj
  #0 to 15 vs. 50 to 70: 0.0682 p < 0.1
  #15 to 30 vs. 50 to 70: 0.0710 p < 0.1
  #0 to 15 vs. 70 to 100: 0.1050
  #15 to 30 vs. 70 to 100: 0.1095
MAOC_change_mmol.let <- multcompLetters4(MAOC_change_mmol.aov, MAOC_change_mmol.tk, 
                                         threshold = 0.1) #adjust threshold
MAOC_change_mmol.let$depth #Fig 3C - 0 to 100: b, b, ab, a, ab (for alpha = 0.1) !!!

MAOM_13C_mmol_plot <- ggplot(zp,aes(x=fct_rev(depth), y=MAOM_13C_mmol))+
  geom_boxplot() +
  theme_zp() + 
  scale_fill_manual(values = zpPalette) +
  labs(y="13C Glucose in MAOC (mmol C)",x="Soil Depth (cm)") +
  scale_y_continuous(limits = c(0,.3)) +
  coord_flip()
MAOM_13C_mmol_plot

MAOM_13Cprop_plot <- ggplot(zp,aes(x=fct_rev(depth), y=MAOM_13Cprop*100))+
  geom_boxplot() +
  theme_zp() + 
  scale_fill_manual(values = zpPalette) +
  labs(y=" MAOC Derived from 13C glucose (%)",x="Soil Depth (cm)") +
  scale_y_continuous(limits = c(0,4)) +
  coord_flip()
MAOM_13Cprop_plot

MAOC_change_mmol_plot <- ggplot(zp,aes(x=fct_rev(depth), y=MAOC_change_mmol))+
  geom_boxplot() +
  theme_zp() + 
  scale_fill_manual(values = zpPalette) +
  labs(y="Change in MAOC with Glucose (mmol C)",x="Soil Depth (cm)") +
  scale_y_continuous(limits = c(-12,2)) +
  coord_flip() +
  geom_hline(yintercept=0, linetype = "dashed")
MAOC_change_mmol_plot

#### Combined SOC Plots - Figure 1 ####

fig_1 <- plot_grid( 
  gC_cm3soilplot + theme(legend.position = "none", text = element_text(size=14),
                         axis.title.y = element_text(size=14)) +
    scale_y_continuous(expand = expansion(mult = c(0,.1))), 
  MAOC_cm3_plot + theme(legend.position = "none", text = element_text(size=14)) +
    labs(y="Carbon in MAOM (g·cm3)",x=NULL) +
    scale_y_continuous(expand = expansion(mult = c(0,.1))), 
  POC_cm3_plot + theme(legend.position = "none", text = element_text(size=14),
                       axis.title.y = element_text(size=14)) +
    scale_y_continuous(expand = expansion(mult = c(0,.1))), 
  MAOC_POCplot + theme(legend.position = "none", text = element_text(size=14)) +
    labs(y="MAOM-C:POM-C",x=NULL) +
    scale_y_continuous(expand = expansion(mult = c(0,.1))),
  nrow = 2,
  ncol = 2,
  align="hv"
) #package cowplot
fig_1

ggsave("Figure_1.png", plot = fig_1, dpi = 300, height = 6,
       width = 7, units = "in")


#### Combined Incubation Plots - Figure 3 ####

fig_3 <- plot_grid(
  MAOM_13C_mmol_plot, 
  MAOM_13Cprop_plot, 
  MAOC_change_mmol_plot,
  nrow = 3,
  ncol = 1,
  align="hv"
) #package cowplot

fig_3

ggsave("Figure_3.png", plot = fig_3, dpi = 300, height = 6,
       width = 3, units = "in")

#### Linear Mixed Effects Models ####

options(na.action = "na.fail") #Set NA option to fail if present

#Identify best model that explains total soil C regardless of horizon

soilClme <- lme(gC_cm3soil ~ AP + BG + NAG + Nit_ug_gsoil_day +
                  Min_ug_gsoil_day + froot_g_cm3soil + R_mmol_nosub, 
                random = ~1|pit/depth,
                data = zp)

#Examine all possible models for soil C
soilClme_output <- dredge(soilClme, 
                          extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))


soilClme_best <- get.models(soilClme_output, 1)[[1]]
soilClme_best


#Model summary for soil C

check_collinearity(soilClme_best) #VIF < 2.5

r.squaredGLMM(soilClme_best) #marginal = 0.72, conditional = 0.98
summary(soilClme_best)
fig2p1 <- as.data.frame(std.coef(soilClme_best, partial = TRUE))[-1,-3]
fig2p1 <- cbind(Fixed_Effect = 
                  c("Fine Root Biomass","NAG Activity","Microbial Respiration"), 
                fig2p1)

#Identify best model that explains total MAOC regardless of horizon

MAOClme <- lme(gMAOC_cm3soil ~ AP + BG + NAG + Nit_ug_gsoil_day +
                  Min_ug_gsoil_day + froot_g_cm3soil + R_mmol_nosub, 
                random = ~1|pit/depth,
                data = zp)

#Examine all possible models for MAOC
MAOClme_output <- dredge(MAOClme, 
                          extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))


MAOClme_best <- get.models(MAOClme_output, 1)[[1]]

#Model summary for soil C

check_collinearity(MAOClme_best) #VIF < 2.5

r.squaredGLMM(MAOClme_best) #marginal = 0.38, conditional = 0.99.. ~1.00
summary(MAOClme_best)
fig2p2 <- as.data.frame(std.coef(MAOClme_best, partial = TRUE))[-1,-3]
fig2p2 <- cbind(Fixed_Effect = 
                  c("Fine Root Biomass","NAG Activity","Microbial Respiration"), 
                fig2p2)

#Identify best model that explains total POC regardless of horizon

POClme <- lme(gPOC_cm3soil ~ AP + BG + NAG + Nit_ug_gsoil_day +
                 Min_ug_gsoil_day + froot_g_cm3soil + R_mmol_nosub, 
               random = ~1|pit/depth,
               data = zp)

#Examine all possible models for POC
POClme_output <- dredge(POClme, 
                         extra = c("R^2", F = function(x) summary(x)$fstatistic[[1]]))


POClme_best <- get.models(POClme_output, 1)[[1]]

#Model summary for POC

r.squaredGLMM(POClme_best) #marginal = 0.39, conditional = 0.99.. ~1.00
summary(POClme_best)

fig2p3 <- as.data.frame(std.coef(POClme_best, partial = TRUE))[-1,-3]
fig2p3 <- cbind(Fixed_Effect = 
                  c("Fine Root Biomass"), 
                fig2p3)

#Figure 2: Std coefficients & St Errors for soil C, MAOC, and POC

fig2.data <- bind_rows(fig2p1, fig2p2)
rm(fig2p1, fig2p2)

fig2.data <- cbind(Variable = 
                  c("Soil C","Soil C","Soil C",
                    "MAOC","MAOC","MAOC"), 
                  fig2.data)

colnames(fig2.data) <- c("Variable", "Fixed_Effect", "Std_Coef", "SE")
str(fig2.data)
fig2.data <- fig2.data %>% 
  mutate(across(c(Variable, Fixed_Effect), as.factor))

fig2.data$Fixed_Effect <- ordered(fig2.data$Fixed_Effect, 
                                  levels = c("Fine Root Biomass",
                                             "NAG Activity",
                                             "Microbial Respiration"))
fig2.data$Variable <- ordered(fig2.data$Variable, 
                                  levels = c("Soil C", "MAOC"))

fig2.bars <- ggplot(fig2.data,aes(x=Fixed_Effect, y=Std_Coef, fill=Variable))+
  geom_bar(position=position_dodge(),stat="identity") +
  geom_errorbar(aes(ymin=Std_Coef-SE,ymax=Std_Coef+SE), 
                width=0.2,position=position_dodge(0.9)) +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.9),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) + 
  scale_fill_manual(values = zpPalette2) + 
  labs(y="Standardized Coefficient", x = element_blank())
fig2.bars

fig2tbl <- as.data.frame(matrix(,2,3))
colnames(fig2tbl) <- c("Variable", "Total Soil C", "MAOC")
fig2tbl$Variable <- c("AIC", "Marginal R2")

fig2tbl[1,2] <- summary(soilClme_best)[[22]]
fig2tbl[1,3] <- summary(MAOClme_best)[[22]]

fig2tbl[2,2] <- r.squaredGLMM(soilClme_best)[1]
fig2tbl[2,3] <- r.squaredGLMM(MAOClme_best)[1]

fig2tbl.plot <- tableGrob(fig2tbl, theme = ttheme_minimal(), rows = NULL)
grid.newpage()
grid.draw(fig2tbl.plot)

fig_2 <- plot_grid(
  fig2.bars, 
  fig2tbl.plot,
  nrow = 2,
  ncol = 1,
  rel_heights = c(4,1)
)

ggsave("fig_2.png", plot = fig_2, dpi = 300, height = 5, width = 5) #also as svg

#Table 2

write.csv(fig2.data, file = "Table_2_StdCoef.csv", row.names = F) #Std coef
#also added the one for POC (from fig2p3) for the full table
