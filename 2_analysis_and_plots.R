## Code accompanying Schmidt et al. (2023) Genetic diversity and IUCN Red List status. Conservation Biology
## 2. Analyses and plots

# Libraries ------
library(tidyr)
library(dplyr)
library(MASS)
library(caret)

# plots
library(patchwork)
library(viridis)
library(extrafont)

# Data ------------
# Mitochondrial DNA (mtDNA)
## Data from: Canteri et al. 2021 Ecography, IUCN Red List protects avian genetic diversity
## https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.05895
mtdna <- read.csv('mtdna_data.csv', header = TRUE) %>% 
  mutate(IUCN = factor(IUCN, levels = c('CR',
                                        'EN',
                                        'VU',
                                        'NT',
                                        'LC')),
         IUCN_fact = as.factor(IUCN_fact))

# MacroPopGen (microsatellites)
## Data from: Lawrence et al. 2019 Scientific Data, Geo-referenced population-specific microsatellite data across American continents, the MacroPopGen Database
## https://www.nature.com/articles/s41597-019-0024-7
mpg <- read.csv('microsat_data.txt', header = TRUE, 
                      sep = '\t')

# Put categories in the correct order:
mpg <- mpg %>% 
  mutate(rlcat = factor(rlcat, levels = c('CR',
                                          'EN',
                                          'VU',
                                          'NT',
                                          'LC',
                                          'DD')))
# Whole genome (WGS) data
## Data from: Brüniche-Olsen et al. 2021 Proceedings B, Life-history traits and habitat availability shape genomic diversity in birds: implications for conservation
## https://royalsocietypublishing.org/doi/10.1098/rspb.2021.1441
wgs <- read.csv("wgs_data.csv", h=T)
wgs$IUCN_fact <- as.factor(wgs$IUCN_fact)
wgs <- wgs %>% 
  mutate(category = factor(category, levels = c('CR',
                                                'EN',
                                                'VU',
                                                'NT',
                                                'LC')))
# mtDNA ---------------
## ordinal model ####
ord_mod <- polr(IUCN_fact ~ GD, data = mtdna, Hess = T)
summary(ord_mod)
predictIUCN <- predict(ord_mod, mtdna, predict = response)
confusionMatrix(mtdna$IUCN_fact, predictIUCN)

## Logistic regression with binary threatened vs non-threatened categories ####
log_mod <- glm(IUCN_bin ~ GD, family = binomial(link=logit), data = mtdna)
summary(log_mod)

# confusion matrix
# predicted probabilities
Yhat <- fitted(log_mod)

# choose a threshold for dichotomizing according to predicted probability
thresh  <- 0.5
YhatFac <- cut(Yhat, breaks=c(-Inf, thresh, Inf), labels=c("0", "1"))
confusionMatrix(as.factor(mtdna$IUCN_bin), YhatFac)


# Macropopgen ----------------

## IUCN cat ####
mpg$IUCN_fact <- as.factor(mpg$IUCN_fact)
mpg$rlcat <- as.factor(mpg$rlcat)
#mpg$IUCN_binary <- as.factor(mpg$IUCN_binary)

# Class data subsets
# Ignore class CEPHALASPIDOMORPHI; only 1 species
mammal <- filter(mpg, class == 'MAMMALIA')
bird <- filter(mpg, class == 'AVES')
reptile <- filter(mpg, class == 'REPTILIA')
amph <- filter(mpg, class == 'AMPHIBIA')
fish <- filter(mpg, class == 'ACTINOPTERYGII')

## ordinal models ####
ord_mam <- polr(IUCN_fact ~ meanGD, data = mammal, Hess = T)
ord_bir <- polr(IUCN_fact ~ meanGD, data = bird, Hess = T)
ord_rep <- polr(IUCN_fact ~ meanGD, data = reptile, Hess = T)
ord_amp <- polr(IUCN_fact ~ meanGD, data = amph, Hess = T)
ord_fis <- polr(IUCN_fact ~ meanGD, data = fish, Hess = T)

summary(ord_mam)
summary(ord_bir)
summary(ord_rep)
summary(ord_amp)
summary(ord_fis)

### confusion matrix ####
predictIUCNm <- predict(ord_bir, mammal, predict = response)
predictIUCNb <- predict(ord_mam, bird, predict = response)
predictIUCNr <- predict(ord_rep, reptile, predict = response)
predictIUCNa <- predict(ord_amp, amph, predict = response)
predictIUCNf <- predict(ord_fis, fish, predict = response)

#confusion matrix
confusionMatrix(bird$IUCN_fact, predictIUCNb)
confusionMatrix(mammal$IUCN_fact, predictIUCNm)
confusionMatrix(reptile$IUCN_fact, predictIUCNr)
confusionMatrix(amph$IUCN_fact, predictIUCNa)
confusionMatrix(fish$IUCN_fact, predictIUCNf)

## logistic regressions ####
logr_bir <- glm(IUCN_binary ~ meanGD, data = bird, family = binomial(link=logit))
logr_mam <- glm(IUCN_binary ~ meanGD, data = mammal, family = binomial(link=logit))
logr_rep <- glm(IUCN_binary ~ meanGD, data = reptile, family = binomial(link=logit))
logr_amp <- glm(IUCN_binary ~ meanGD, data = amph, family = binomial(link=logit))
logr_fis <- glm(IUCN_binary ~ meanGD, data = fish, family = binomial(link=logit))

Yhat <- fitted(logr_bir)
predictIUCNbl <- fitted(logr_bir)
YhatFac <- cut(predictIUCNbl, breaks=c(-Inf, thresh, Inf), labels=c("0", "1"))
confusionMatrix(as.factor(bird$IUCN_binary), YhatFac)

Yhat <- fitted(logr_mam)
predictIUCNml <- fitted(logr_mam)
YhatFac <- cut(Yhat, breaks=c(-Inf, thresh, Inf), labels=c("0", "1"))
confusionMatrix(as.factor(mammal$IUCN_binary), YhatFac)

Yhat <- fitted(logr_amp)
predictIUCNal <- fitted(logr_amp)
YhatFac <- cut(Yhat, breaks=c(-Inf, thresh, Inf), labels=c("0", "1"))
confusionMatrix(as.factor(amph$IUCN_binary), YhatFac)

Yhat <- fitted(logr_rep)
predictIUCNrl <- fitted(logr_rep)
YhatFac <- cut(Yhat, breaks=c(-Inf, thresh, Inf), labels=c("0", "1"))
confusionMatrix(as.factor(reptile$IUCN_binary), YhatFac)

Yhat <- fitted(logr_fis)
predictIUCNfl <- fitted(logr_fis)
YhatFac <- cut(Yhat, breaks=c(-Inf, thresh, Inf), labels=c("0", "1"))
confusionMatrix(as.factor(fish$IUCN_binary), YhatFac)

# WGS ####
## ordinal model ####
ord_modw <- polr(IUCN_fact ~ He, data = wgs, Hess = T)
summary(ord_modw)
predictIUCNw <- predict(ord_modw, wgs, predict = response)
confusionMatrix(wgs$IUCN_fact, predictIUCNw)

## Logistic regression with binary threatened vs non-threatened categories ####
log_modw <- glm(IUCN_binary ~ He, family = binomial(link=logit), data = wgs)
summary(log_modw)

# confusion matrix
# predicted probabilities
Yhat <- fitted(log_modw)
predictIUCNwl <- fitted(log_modw)
# choose a threshold for dichotomizing according to predicted probability
thresh  <- 0.5
YhatFac <- cut(Yhat, breaks=c(-Inf, thresh, Inf), labels=c("0", "1"))
confusionMatrix(as.factor(wgs$IUCN_binary), YhatFac)


# Plots ----------------------
## boxplots ####
######### fish,     amphibians, birds,    lamprey,   mammals,   reptiles
pal <- c('#ED8B16', '#E1523D', '#C2BB00', '#734100', '#003547', '#0B8301')

class_box <- ggplot(data = mpg %>% filter(class != 'CEPHALASPIDOMORPHI'), aes(y = meanGD, x = rlcat)) + 
  geom_boxplot(notch = FALSE, lwd = 0.8, outlier.shape = NA) +
  geom_jitter(aes(fill = class), width = 0.08, alpha = 0.7, shape = 21) +
  scale_fill_manual(values = pal,
                    guide = 'none') +
  labs(x = '', y = 'gene diversity') +
  facet_wrap(~class) +
  theme_minimal() + 
  theme(text=element_text(size=16, 
                          family="Roboto"))

mt_box <- ggplot(data = mtdna, aes(y = GD, x = IUCN)) + 
  geom_boxplot(notch = FALSE, lwd = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.08, fill = '#C2BB00', alpha = 0.7, shape = 21) +
  labs(x = '', y = 'nucleotide diversity') +
  theme_minimal() + 
  theme(text=element_text(size=16, 
                          family="Roboto"))

wgs_box <- ggplot(data = wgs, aes(y = He, x = category)) + 
  geom_boxplot(notch = FALSE, lwd = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.08, fill = '#C2BB00', alpha = 0.7, shape = 21) +
  labs(x = '', y = 'observed heterozygosity') +
  theme_minimal() + 
  theme(text=element_text(size=16, 
                          family="Roboto"))

layout <- "
AABB
CCCC
CCCC
"

(f3<-mt_box + wgs_box + class_box + plot_annotation(tag_levels = 'A') + plot_layout(design = layout))

#ggsave(f3, filename = "fig3.pdf", device = cairo_pdf, width = 10, height = 8, units = "in")


## mtdna vs microsat div ####
mtdna$SP <- gsub("_", " ", mtdna$SP)
share_sp <- which(mpg$scientific_name %in% c(mtdna$SP))

share_mpg <- mpg[share_sp,]
share_mtd <- mtdna %>%
  filter(SP %in% c(share_mpg$scientific_name))

micro_mtd <- merge(share_mpg[,c('scientific_name', 'meanGD', 'rlcat')],
                   share_mtd[,c('SP', 'GD')],
                   by.x = 'scientific_name', by.y = 'SP')

micro_mtd$rlcat <- factor(micro_mtd$rlcat, levels = c('LC', 'NT', 'VU', 'EN', 'CR'))
micro_mtd$rowID <- seq(1:nrow(micro_mtd))

(f2 <- ggplot(micro_mtd, aes(y = GD, x = meanGD, color = rlcat)) + 
  geom_point(size = 2) +
  scale_color_viridis(discrete = T, option = "inferno", 
                      begin = 0.15, end = 0.8, name = NULL) +
  labs(y = "mitochondrial genetic diversity",
       x = "microsatellite genetic diversity") +
  theme_minimal() +
  theme(text=element_text(size=12, 
                          family="Roboto")))
ggsave(f2, filename = "fig2.pdf", device = cairo_pdf, width = 6, height = 5, units = "in")

cor.test(micro_mtd$meanGD, micro_mtd$GD)



# Confusion matrices ---------
confused_plot <- function(predicted, observed, title){
  if(class(predicted) == 'factor'){
    hmcm <- as.data.frame(table(predicted, observed)/sum(table(observed)))
  }
  else{
    YhatFac <- cut(predicted, breaks=c(-Inf, 0.5, Inf), labels=c("0", "1"))
    hmcm <- as.data.frame(table(YhatFac, observed)/sum(table(observed)))
  }
  names(hmcm) <- c('predicted', 'observed', 'Freq')
  hmcm$perc <- round(hmcm$Freq*100, 2)
  
  ggplot(data = hmcm, aes(y = predicted, x = observed, fill = perc, label = perc)) +
    geom_tile(color = '#FFFFFF') +
    geom_text(color = '#000000', family = 'Roboto Black') +
    scale_fill_gradientn(colours = c('#FFFFFF', '#c296b7ff'), 
                         limits=c(0, 100), name = "Frequency (%)") +
    scale_y_discrete(limits=rev) +
    labs(title = title) +
    theme_minimal() +
    theme(text=element_text(size=12, 
                            family="Roboto"))
}

## ordinal models -----
mtDNA_ordinal <- confused_plot(predictIUCN, mtdna$IUCN_fact, "mtDNA")
mtDNA_log <- confused_plot(predict, mtdna$IUCN_bin, "mtDNA")
mam_ordinal <- confused_plot(predictIUCNm, mammal$IUCN_fact, "microsatellites (mammals)")
mam_log <- confused_plot(predictIUCNml, mammal$IUCN_binary, "microsatellites (mammals)")
bir_ordinal <- confused_plot(predictIUCNb, bird$IUCN_fact, "microsatellites (birds)")
bir_log <- confused_plot(predictIUCNbl, bird$IUCN_binary, "microsatellites (birds)")
rep_ordinal <- confused_plot(predictIUCNr, reptile$IUCN_fact, "microsatellites (reptiles)")
rep_log <- confused_plot(predictIUCNrl, reptile$IUCN_binary, "microsatellites (reptiles)")
amp_ordinal <- confused_plot(predictIUCNa, amph$IUCN_fact, "microsatellites (amphibians)")
amp_log <- confused_plot(predictIUCNal, amph$IUCN_binary, "microsatellites (amphibians)")
fis_ordinal <- confused_plot(predictIUCNf, fish$IUCN_fact, "microsatellites (fish)")
fis_log <- confused_plot(predictIUCNfl, fish$IUCN_binary, "microsatellites (fish)")
wgs_ordinal <- confused_plot(predictIUCNw, wgs$IUCN_fact, "WGS")
wgs_log <- confused_plot(predictIUCNwl, wgs$IUCN_binary, "WGS")


# mtDNA & WGS
#png('confus_WGS_mt.png', width = 8, height = 8, res = 300, units = 'in')
(mtDNA_ordinal + mtDNA_log)/(wgs_ordinal + wgs_log) + plot_layout(guides = 'collect')
#dev.off()

# microsats
#png('confus_msat.png', width = 8, height = 10, res = 300, units = 'in')
(fis_ordinal + fis_log)/(amp_ordinal + amp_log)/
  (bir_ordinal + bir_log)/(mam_ordinal + mam_log)/
  (rep_ordinal + rep_log) + plot_layout(guides = 'collect')
#dev.off()
