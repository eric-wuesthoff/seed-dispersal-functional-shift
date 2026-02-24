###Maromizaha seed disperal functional shift project
###This code represents analysis for adult plant communities (part 1),
###lemur communities (part 2), seed rain (part 3), dispersed seeds (part 4),
###seed size (part 5)

###Manuscript authors:Eric F Wuesthoff, Lovasoa Natenaina Ratolojanahary, 
###M. Rahelison, Tojoharilala Fenohasina Merison, Harison Rabarison, Amy E. Dunham

#PREFACE Load packages####

library(ggplot2)
library(ggpubr)
library(ggeffects)
library(ggbreak)
library(TMB)
library(vegan)
library(glmmTMB)
library(cplm)
library(Distance)
library(dplyr)
library(dtplyr)
library(tidyr)
library(glmmTMB)
library(lme4)
library(generics)
library(ordinalNet)
library(ggsignif)
library(ggpubr)
library(grid)
library(gridExtra)
library(emmeans)

#-------------------------------------------------------------------------------------------------------------------------------------------------------
#PART 1: Plant Community Diversity & Fruiting####

#1a: Load & configure data####

#taxa abundance by trap cluster
overhang <- read.csv ('Final Data/plant_community_overhang_mmz2022.csv')
#combine trap clusters for one row per transect
overhang_sum <- rowsum(x = overhang[c(1:43),c(7:181)], group = overhang$Transect[c(1:43)])
Transect <- c("a","b","c","d","e","f")
Habitat <- c("protected","protected","savoka","restoration","restoration","protected")
Regeneration <- c("intact","intact","regenerating","regenerating","regenerating","intact")
Anthro.rank <- c("2","1","6","5","4","3")
overhang_sum <- data.frame(Transect, Habitat, Regeneration, Anthro.rank, overhang_sum)

#plant flowering and fruiting rates by survey number
overhang.phenology <- read.csv('Final Data/overhang_phenology.csv')
#create vector of proportion of plants that are non-native (1 sig fig)
proportion.nonnative <-c(0.004, 0.004, 0.3, 0.4, 0.3, 0.05)

#2a Calculate diversity metrics####

#using renyi() function from vegan package
overhang_ren <- renyi(overhang_sum[5:179], scale = c(0,1,2,Inf), hill = TRUE)
colnames(overhang_ren) <- (c("richness", "shannons","inv_simpson","weighted_div")) 

overhang_ren <- merge(overhang_ren, overhang_sum, by = 0) %>%
  mutate(overhang_ren, Transect = Row.names) 
overhang_ren$total_abundance <-apply(overhang_ren[10:184], 1, sum)

overhang_ren <- data.frame(overhang_ren,proportion.nonnative)

#testing abundance/diversity metrics for normality
shapiro.test(overhang_ren $total_abundance) #p=0.01759, total abundance non-normal
hist(overhang_ren $total_abundance)

shapiro.test(overhang_ren $richness) #p=0.04147, richness non-normal
hist(overhang_ren $richness)

shapiro.test(overhang_ren $shannons) #p=0.04564, shannon's diversity non-normal
hist(overhang_ren $shannons)

shapiro.test(overhang_ren $inv_simpson) #p=0.05723, inverse simpson's non-normal
hist(overhang_ren $inv_simpson)

shapiro.test(overhang_ren$proportion.nonnative) #p=0.1072, proportion non-native plants normal
hist(overhang_ren$proportion.nonnative)

#3c: Model differences in plant community diversity between mature & regenerating forests####

#glm used (as opposed to glmer) due to singular boundary fit using GLMM
#Models do not hold transect as random effect
#Summarized in Table S1 in article supplement
summary(glm(richness~Regeneration,overhang_ren,family=poisson(link="log"))) #significantly lower in regen
summary(glm(shannons~Regeneration,overhang_ren,family=Gamma(link="log"))) #significantly lower in regen
summary(glm(inv_simpson~Regeneration,overhang_ren,family=Gamma(link="log"))) #significantly lower in regen
summary(lm(proportion.nonnative~Regeneration,overhang_ren)) #significantly higher in regen

#1d: Visualize mean plant community metrics between forest types####

overhang_richness_plot<-
  ggplot(overhang_ren,aes(x=Regeneration, y=richness,fill=Regeneration))+
  geom_boxplot(show.legend=FALSE)+
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  scale_x_discrete(label = c("Mature", "Regenerating"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle=15,vjust=0.6))+
  theme(plot.background = element_rect(fill='transparent',color='transparent'))+
  labs(y = ("Richness"), x = "")+
  annotate(geom = "text", x = 1.5, y = 63, label = "***", cex = 6)+
  annotate(geom = "text", x = 0.5, y = 90, label = "A", cex = 8)

overhang_shannon_plot<-
  ggplot(overhang_ren,aes(x=Regeneration, y=shannons,fill=Regeneration))+
  geom_boxplot(show.legend=FALSE)+
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  scale_x_discrete(label = c("Mature", "Regenerating"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle=15,vjust=0.6))+
  theme(plot.background = element_rect(fill='transparent',color='transparent'))+
  labs(y = ("Shannon Index"), x = "")+
  annotate(geom = "text", x = 1.5, y = 46, label = "***", cex = 6)+
  annotate(geom = "text", x = 0.5, y = 72, label = "B", cex = 8)

overhang_evenness_plot<-
  ggplot(overhang_ren,aes(x=Regeneration, y=inv_simpson,fill=Regeneration))+
  geom_boxplot(show.legend=FALSE)+
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  scale_x_discrete(label = c("Mature", "Regenerating"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle=15,vjust=0.6))+
  theme(plot.background = element_rect(fill='transparent',color='transparent'))+
  labs(y = ("Inverse Simpson Index"), x = "")+
  annotate(geom = "text", x = 1.5, y = 40, label = "***", cex = 6)+
  annotate(geom = "text", x = 0.5, y = 59, label = "C", cex = 8)

overhang_nonnative_plot<-
  ggplot(overhang_ren,aes(x=Regeneration, y=proportion.nonnative,fill=Regeneration))+
  geom_boxplot(show.legend=FALSE)+
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  scale_x_discrete(label = c("Mature", "Regenerating"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle=15,vjust=0.6))+
  theme(plot.background = element_rect(fill='transparent',color='transparent'))+
  labs(y = ("Propotion of\nPlants Nonnative"), x = "")+
  annotate(geom = "text", x = 1.5, y = 0.2, label = "***", cex = 6)+
  annotate(geom = "text", x = 0.47, y = 0.4, label = "D", cex = 8)


#Figure S1 in article supplement
#Combine richness, Shannon, & inverse Simpson figures
grid.arrange(overhang_richness_plot,overhang_shannon_plot, overhang_evenness_plot,overhang_nonnative_plot,ncol=2)


#1e: Model differences in fruiting proportions between forest types & over time####
#Models summarized in Table S1 of article supplement

#All plants
#test for normality
shapiro.test(overhang.phenology$Proportion_indiv_all_fruit_total)
hist(overhang.phenology$Proportion_indiv_all_fruit_total)

#Beta distribution chosen due to positive proportion data
summary(glmmTMB(Proportion_indiv_all_fruit_total~Regeneration+Survey.week+(1|Transect),
                overhang.phenology,
                family=beta_family(link="logit"))) 
#Findings:
#Regenerating forests have significantly higher proportion of fruiting individuals
#Proportion of fruiting plants decreases over time

#Non-native plants
#test for normality
shapiro.test(overhang.phenology$Proportion_indiv_all_fruit_nonnative)
hist(overhang.phenology$Proportion_indiv_all_fruit_nonnative)

#Tweedie distribution chosen due to high zeros in continuous proportion data
summary(glmmTMB(Proportion_indiv_all_fruit_nonnative~Regeneration+Survey.week+(1|Transect),
                overhang.phenology,
                family=tweedie)) 
#Findings:
#Regeneratiing forests have significantly higher proportions of fruiting nonnative plants
#Proportion of fruiting nonnatives does not significantly decrease over time (p=0.0719)

#1f: Visualize plant community fruiting between forest types & over time####
total.fruit.proportion.plot <- ggplot(overhang.phenology,aes(x=factor(Survey.week,ordered=F),
                                                       y=Proportion_indiv_all_fruit_total,fill=Regeneration))+
  geom_smooth(aes(group=Regeneration,color=Regeneration),method="glm",show.legend = FALSE) +
  geom_boxplot(show.legend=FALSE)+
  #geom_point()+
  scale_fill_manual(values = c("slateblue", "skyblue"))+
  scale_color_manual(values = c("slateblue", "skyblue"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=16,face="bold"),
        axis.title.x=element_blank(),
        plot.margin = unit(c(0,0,1.5,1.5), "cm"))+
  ylim(0,0.53)+
  labs(y = ("Proportion of\nTotal Plants Fruiting"))+
  annotate(geom = "text", x = 0.75, y = 0.5, label = "a", cex = 8)


nonnative.fruit.proportion.plot <- ggplot(overhang.phenology,aes(x=factor(Survey.week,ordered=F),
                                                       y=Proportion_indiv_all_fruit_nonnative,fill=Regeneration))+
  #geom_smooth(aes(group=Regeneration,color=Regeneration),method="glm",show.legend = FALSE) +
  geom_boxplot(show.legend=FALSE)+
  #geom_point()+
  scale_fill_manual(values = c("slateblue", "skyblue"))+
  scale_color_manual(values = c("slateblue", "skyblue"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=16,face="bold"),
        axis.title.x=element_blank(),
        plot.margin = unit(c(0,0,1.5,1.5), "cm"))+
  ylim(0,0.53)+
  labs(y = ("Proportion of\nNonnative Plants Fruiting"))+
  annotate(geom = "text", x = 0.75, y = 0.5, label = "b", cex = 8)

#combine total & non-native fruiting plots
#Figure 4 a,b in article
fruit.plots <- grid.arrange(total.fruit.proportion.plot,nonnative.fruit.proportion.plot,ncol=2, 
                            top=textGrob("Fruiting Plants",gp=gpar(fontsize=20,font=2),x=0,hjust=0,vjust=-0.2))


#-------------------------------------------------------------------------------------------------------------------------------------------------------
#PART 2: Lemur Transect Surveys (Encounter Rates)####

#2a: Load & configure data#####
lemur_species_encounters<- read.csv("Final Data/lemur_survey_species_encounter_rates.csv")
#Lemur diet/activity guilds:
#1. diurnal frugivores: Vareica, Eulemur
#2. diurnal folivores: Hapalemur, Propithecus, Indri
#3. nocturnal folivores: Lepilemur, Avahi
#4. nocturnal omnivores: Microcebus

#encounter rates are for each transect are calculated by (number of animals observed) / (length of transect)

lemur_encounter.rates <-read.csv("Final Data/lemur_functional_encounter_rates.csv")
#test for normality
shapiro.test(lemur_encounter.rates$func.encounter.rate) #p-value < 2.2e-16, non-normal

##2b: Model lemur encounter rates by diet/activity guild ####

#Model summarized in Table S3 of article supplment
simp.funcmod1 <-glmmTMB(func.encounter.rate~regeneration+(1|transect), 
                data=lemur_encounter.rates,family=tweedie)
simp.funcmod2 <- glmmTMB(func.encounter.rate~functional.group+(1|transect), 
                data=lemur_encounter.rates,family=tweedie)

funcmod1 <- glmmTMB(func.encounter.rate~functional.group+regeneration+(1|transect), 
                data=lemur_encounter.rates,family=tweedie)
summary(funcmod1)
#Findings: microcebus encountered sig more than diurnal folivores & frugivores

funcmod2 <- glmmTMB(func.encounter.rate~functional.group+regeneration+functional.group*regeneration+(1|transect), 
                    data=lemur_encounter.rates,family=tweedie) #interaction between diurnal frugivores and regenerating habitat
summary(funcmod2)

anova(simp.funcmod1, simp.funcmod2, funcmod1, funcmod2, test="LRT")


#pairwise comparisons
em <- emmeans(funcmod1, "functional.group")
contrast(em, "pairwise", adjust = "Tukey")

##2c Model lemur encounter rates by body size functional group ####

lemur.size <- read.csv("Final Data/lemur_size_combined.csv")
lemur.size$survey.ID<-as.factor(lemur.size$survey.ID)
#test for normality
shapiro.test(lemur.size$size.encounter.rate) #p-value < 2.2e-16, non-normal


simple.sizemod1.reorder <- glmmTMB(size.encounter.rate~regeneration+(1|transect), 
                data=lemur.size,family=tweedie)

simple.sizemod2.reorder <- glmmTMB(size.encounter.rate~body.size.reorder+(1|transect), 
                                           data=lemur.size,family=tweedie)

sizemod1.reorder <-glmmTMB(size.encounter.rate~body.size.reorder+regeneration+(1|transect), 
                           data=lemur.size,family=tweedie)

sizemod2.reorder <-
  glmmTMB(size.encounter.rate~body.size.reorder+regeneration+body.size.reorder*regeneration+(1|transect), 
          data=lemur.size,family=tweedie)

anova(simple.sizemod1.reorder, simple.sizemod2.reorder, sizemod1.reorder, sizemod2.reorder, test="LRT")

summary(sizemod2.reorder)

#2d Model small (mouse lemur) primate encounter rates by forest type ####

noct.lemur.size <-lemur.size[lemur.size$period=="Night",]
micro.size <- noct.lemur.size[noct.lemur.size$body.size=="small",]

shapiro.test(micro.size$size.encounter.rate)
hist(micro.size$size.encounter.rate)

smallmod <-
  lm(size.encounter.rate~regeneration,
          data=micro.size)
summary(smallmod)
aov(size.encounter.rate~regeneration,dat=micro.size)

#2e: Lemur heatmap visualizations ####

#cumulative encounter rate by species
lemur.abun <- read.csv("Final Data/lemur_abundances.csv")
lemur.abun$species <- factor(lemur.abun$species,
                             levels=c("ms","mr","ml","lm","al","hg","pd","ii","ef","vv"))
lemur.abun$anthro.rank <- factor(lemur.abun$anthro.rank,
                                 levels=c("MATURE1","MATURE2","MATURE3","REGEN1","REGEN2","REGEN3"))
lemur.abun$rate <- lemur.abun$rate*100

#plot species heatmap
lemur.sp.abundance<-
  ggplot(lemur.abun,aes(x=species,y=anthro.rank,fill=rate))+
  geom_tile(color="black")+
  geom_text(aes(label=round(rate,digits=1)))+
  coord_flip()+scale_fill_gradient(low="yellow",high="red",limits=c(0,4))+
  #scale_x_discrete(position = "top")+
  theme_classic()+labs(y="Transects",x="",fill="Cumulative Encounter Rate\n(Individuals per 100 Meters)")+
  theme(axis.title.y=element_text(angle=0,vjust=0.5,size=14),legend.position=c(1.27,-0.18),#legend.position = c(0.8,0.45),
        legend.direction = "horizontal",legend.title.position = "top",
        axis.text=element_text(size=12),axis.text.y=element_text(face="italic"), axis.text.x=element_text(angle=25,vjust=0.5),
        legend.title=element_text(size=10,hjust=0.5),legend.text=element_text(size=10), 
        legend.background = element_rect(linetype="solid",colour = "black"),
        axis.title.x=element_text(size=14,margin=margin(t=11,b=0,r=0,l=0)),
        plot.margin = unit(c(0,1,1,0), "cm"))+
  scale_x_discrete(position = "top",labels=c("Microcebus lehilahytsara","Lepilemur mustelinus","Avahi laniger","Hapalemur griseus",
                                             "Propithecus diadema","Indri indri","Eulemur fulvus","Varecia variegata"))

##cumulative encounter rate by diet/activity guild 
lemur.enc.transect2 <- read.csv("Final Data/lemur_transect_rates2.csv")
lemur.enc.transect2$functional.group <- factor(lemur.enc.transect2$functional.group,
                                               levels=c("diu.frug","diu.fol","noc.fol","micro"))
lemur.enc.transect2$transect <- factor(lemur.enc.transect2$transect,
                                       levels=c("m1","m2","m3","r1","r2","r3"))
lemur.enc.transect2 <-lemur.enc.transect2[lemur.enc.transect2$rate!=0,]
lemur.enc.transect2$rate <-lemur.enc.transect2$rate*100

#plot diet/activity guild  heatmap
lemur.func.abundance<-
  ggplot(lemur.enc.transect2,aes(x=functional.group,y=transect,fill=rate))+
  geom_tile(color="black")+
  geom_text(aes(label=round(rate,digits=1)))+
  coord_flip()+scale_fill_gradient(na.value="white",low="yellow",high="red",limits=c(0,4))+
  #scale_x_discrete(position = "top")+
  theme_classic()+labs(y="Transects",x="",fill="Cumulative Encounter Rate\n(Individuals per 100 Meters)")+
  theme(legend.position = 'none',
        axis.title.y=element_text(angle=0,vjust=0.5,size=14),#legend.position = c(0.8,0.5),
        axis.text=element_text(size=12),axis.text.y=element_text(face="italic"), 
        axis.text.x=element_text(angle=25,vjust=0.5),
        #legend.title=element_text(size=9),legend.text=element_text(size=9),
        axis.title.x=element_text(size=14,margin=margin(t=11,b=0,r=0,l=0)),
        plot.margin = unit(c(0,0,1,0), "cm"))+
  scale_y_discrete(labels=c("MATURE1","MATURE2","MATURE3","REGEN1","REGEN2","REGEN3"))+
  scale_x_discrete(position="top",
                   labels=c("Nocturnal omniovre\n(Microcebus)","Nocturnal folivore","Dirunal folivore",
                            "Dirunal frugivore"),
                   limits=c("micro","diu.fol","noc.fol","diu.frug"))

#Figure 2 in article
#combined plot
grid.arrange(lemur.sp.abundance,lemur.func.abundance,ncol=2) #1500 x 400

#2f: Microcebus encounter boxplot visualization ####

microcebus_encounter <- lemur_encounter.rates %>% filter(period =="Night", functional.group=="_microcebus")

microcebus_encounter$order.transect <- factor(microcebus_encounter$order.transect,
                                 levels=c("MATURE1","MATURE2","MATURE3","REGEN1","REGEN2","REGEN3",order=T))

#Figure 3 in article
ggplot(microcebus_encounter,aes(x=order.transect,y=size.encounter.rate*100,fill=regeneration))+ #microcebus encounters plot
  geom_boxplot()+
  #geom_point(color="purple")
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  scale_x_discrete(label = c("MATURE1", "MATURE2","MATURE3","REGEN1","REGEN2","REGEN3"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle=15,vjust=0.6),
        plot.background = element_rect(fill='transparent',color='transparent'),
        legend.title=element_text(size=12,hjust=0.5),legend.text=element_text(size=12),
        legend.position=c(0.25,0.65),
        legend.background = element_rect(linetype="solid",colour = "black"))+
  labs(y = ("Mouse Lemur Encounter Rate\n(Indivduals per 100 meters)"), x = "Transect",fill="Forest Type")

#-------------------------------------------------------------------------------------------------------------------------------------------------------

#PART 3: Total Seed Rain Diversity & Abundance####

#3a: Load & configure data####

#counts of taxon-identified seeds from seed rain by each transect
rain_counts_species <-read.csv('Final Data/seed_rain_counts_species.csv')
rain_counts_species$transect <- factor(rain_counts_species$transect)
#data on individual seeds recorded in seed rain
rain_individuals <- read.csv('/Users/ericwuesthoff/Desktop/Maromizaha Dispersal (Ch 2)/Spreadsheets/Formatted for analysis/rain.seed.individuals.csv')
#data seed rain & dispersal rates on each transect by week (and associated fruiting rates)
rain.dispersal_rates <-read.csv('Final Data/seed.rain.dispersal.rates.csv')

#3b: Calculate diversity metrics####

#using renyi() function from vegan package
seed.rain_ren <- renyi(rain_counts_species[6:25], scale = c(0,1,2,Inf), hill = TRUE)
colnames(seed.rain_ren) <- (c("richness", "shannons","inv_simpson","weighted_div")) 
seed.rain_ren <- merge(seed.rain_ren, rain_counts_species, by = 0) 

#testing for normality
shapiro.test(seed.rain_ren$richness) #p=0.121, normal
hist(seed.rain_ren$richness)
shapiro.test(seed.rain_ren$shannons) #p=0.361, normal
hist(seed.rain_ren$shannons)
shapiro.test(seed.rain_ren$inv_simpson) #p=0.528, normal

#3c: Model differences in seed rain diversity between mature & regenerating forests ####
#Models summarized in Table S1 of article supplement
summary(lm(richness~regeneration,seed.rain_ren)) #p=0.759
summary(lm(shannons~regeneration,seed.rain_ren)) #p=0.165
summary(lm(inv_simpson~regeneration,seed.rain_ren)) #p=0.087
#Findings:
#Seed rain diversity does not significantly between mature and regenerating forests.

#3d: Visualize mean seed rain community metrics between forest types ####
#richness
seed.rain_richness_plot<-
  ggplot(seed.rain_ren,aes(x=regeneration, y=richness,fill=regeneration))+
  geom_boxplot(show.legend=FALSE)+
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  scale_x_discrete(label = c("Mature", "Regenerating"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle=15,vjust=0.6))+
  theme(plot.background = element_rect(fill='transparent',color='transparent'))+
  labs(y = ("Richness"), x = "")+
  annotation_custom(grid::textGrob(x = unit(0.5,"npc"), y = unit(0.9,"npc"), label = "NS", gp= gpar(cex=1.5)))+
  annotation_custom(grid::textGrob(x = unit(0.05,"npc"), y = unit(0.93,"npc"), label = "a", gp= gpar(cex=2)))

#Hill-Shannon diversity
seed.rain_shannon_plot<-
  ggplot(seed.rain_ren,aes(x=regeneration, y=shannons,fill=regeneration))+
  geom_boxplot(show.legend=FALSE)+
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  scale_x_discrete(label = c("Mature", "Regenerating"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle=15,vjust=0.6))+
  theme(plot.background = element_rect(fill='transparent',color='transparent'))+
  labs(y = ("Shannon Index"), x = "")+
  annotation_custom(grid::textGrob(x = unit(0.5,"npc"), y = unit(0.9,"npc"), label = "NS", gp= gpar(cex=1.5)))+
  annotation_custom(grid::textGrob(x = unit(0.05,"npc"), y = unit(0.93,"npc"), label = "b", gp= gpar(cex=2)))

#inverse Simpson diversity
seed.rain_evenness_plot<-
  ggplot(seed.rain_ren,aes(x=regeneration, y=inv_simpson,fill=regeneration))+
  geom_boxplot(show.legend=FALSE)+
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  scale_x_discrete(label = c("Mature", "Regenerating"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle=15,vjust=0.6))+
  theme(plot.background = element_rect(fill='transparent',color='transparent'))+
  labs(y = ("Inverse Simpson Index"), x = "")+
  annotation_custom(grid::textGrob(x = unit(0.5,"npc"), y = unit(0.9,"npc"), label = "NS", gp= gpar(cex=1.5)))+
  annotation_custom(grid::textGrob(x = unit(0.05,"npc"), y = unit(0.93,"npc"), label = "c", gp= gpar(cex=2)))

#Figure S3 a,b,c in article supplement
rain.diversity<- grid.arrange(seed.rain_richness_plot,seed.rain_shannon_plot, seed.rain_evenness_plot,
                                   ncol=3,
                              top=textGrob("Seed Rain",gp=gpar(fontsize=20,font=2)))


#3e: Model differences in seed rain abundance between forest types & over time####

#testing for normality
shapiro.test(rain.dispersal_rates$rain.encounter.rate) #p=1.395e-11, rates non-normal
hist(rain.dispersal_rates$rain.encounter.rate)

shapiro.test(rain.dispersal_rates$nonnative.rain.encounter.rate) #p=1.834e-12, rates non-normal
hist(rain.dispersal_rates$nonnative.rain.encounter.rate)

#GLMMs chosen to account generalize for non-normal data and set transect as random effect.
#Tweedie distribution selected for positive continuous proportions with higher number of zeros

#Total seed rain models

# total rain rates without fruiting
summary(glmmTMB(rain.encounter.rate~regeneration+survey.number +(1|transect),
                family=tweedie,data=rain.dispersal_rates)) #forest type & week significant

# total rain rates with fruiting
summary(glmmTMB(rain.encounter.rate~regeneration+survey.number +fruiting.proportion+(1|transect),
                family=tweedie,data=rain.dispersal_rates)) #fruiting & week significant (forest type NOT sig)

#Findings:
#Seed rain decreases over time during our 7 week study
#Regenerating forests have increased seed rain compared to mature forests
#Effects of forest type are explained by higher proportions fo fruiting plants in regenerating forests

#Non-native seed rain models
#Models summarized in Table 2 of article and Table S2 in supplement

# non-native rain rates without fruiting
summary(glmmTMB(nonnative.rain.encounter.rate~regeneration+survey.number +(1|transect),
                family=tweedie,data=rain.dispersal_rates)) #forest type & week significant

# non-native rain rates with fruiting
summary(glmmTMB(nonnative.rain.encounter.rate~regeneration+survey.number + nonnative.fruiting.proportion+(1|transect),
                family=tweedie,data=rain.dispersal_rates)) #fruiting, week, & forest type significant

#Findings:
#Non-native seed rain abundance decreases over time during our 7 week study.
##Regenerating forests have increased seed rain compared to mature forests
#Effect of forest type remains significant even when controlling for non-native fruiting

#3f: Visualize seed rain abundance between forest types & over time ####

total_rain_plot <- 
  ggplot(rain.dispersal_rates,aes(x=as.factor(survey.number), rain.encounter.rate,fill=regeneration))+
  # geom_point()+
  geom_smooth(aes(group=regeneration,color=regeneration),method="glm",show.legend = FALSE)+
  geom_boxplot(show.legend=FALSE)+ 
  scale_fill_manual(values = c("slateblue", "skyblue"))+
  scale_color_manual(values = c("slateblue", "skyblue"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=16,face="bold"),
        plot.margin = unit(c(0,0,1.5,1.5), "cm"),
        axis.title.x=element_blank(),
  )+ ylim(0,18)+
  labs(y = ("Total Seed Rain\n(Seeds per Sampling Effort)"))+
  #geom_smooth(aes(group=Regeneration,color=Regeneration),method="lm",show.legend = FALSE)
  annotate(geom = "text", x = 0.75, y=17, label = "c", cex = 8)


non.native_rain_plot <- 
  ggplot(rain.dispersal_rates,aes(x=as.factor(survey.number), y=nonnative.rain.encounter.rate,fill=regeneration))+
  geom_smooth(aes(group=regeneration,color=regeneration),method="glm",show.legend = FALSE)+
  geom_boxplot(show.legend=FALSE)+
  scale_fill_manual(values = c("slateblue", "skyblue"))+
  scale_color_manual(values = c("slateblue", "skyblue"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=16,face="bold"),
        plot.margin = unit(c(0,0,1.5,1.5), "cm"),
        axis.title.x=element_blank(),
  )+ ylim(0,18)+
  labs(y = ("Non-native Seed Rain\n(Seeds per Sampling Effort)"))+
  #geom_smooth(aes(group=Regeneration,color=Regeneration),method="lm",show.legend = FALSE)
  annotate(geom = "text", x = 0.75, y =17, label = "d", cex = 8)

#Figure 4 c,d in article
rain.plots <- grid.arrange(total_rain_plot,non.native_rain_plot,ncol=2,
             top=textGrob("Seed Rain",gp=gpar(fontsize=20,font=2),x=0,hjust=0,vjust=-0.2))

#-------------------------------------------------------------------------------------------------------------------------------------------------------
#PART 4: Dispersed Seeds Diversity & Abundance####

#4a: load & configure data####

#counts of taxon-identified dispersed seeds by each transect
dispersal_counts_species <-read.csv ('Final Data/seed_dispersal_counts_species.csv')
#data on individual seeds recorded as dispersed
dispersal_individuals <- read.csv('Final Data/dispersal.seed.individuals.csv')

#4b: Calculate diversity metrics####

#using renyi() function from vegan package
dispersal_ren <- renyi(dispersal_counts_species[6:25], scale = c(0,1,2,Inf), hill = TRUE)
colnames(dispersal_ren) <- (c("richness", "shannons","inv_simpson","weighted_div")) 
dispersal_ren <- merge(dispersal_ren, dispersal_counts_species, by = 0) 

#test for normality
shapiro.test(dispersal_ren$richness) #p=0.70, normal
hist(dispersal_ren$richness)
shapiro.test(dispersal_ren$shannons) #p=0.484, normal
hist(dispersal_ren$shannons)
shapiro.test(dispersal_ren$inv_simpson) #p=0.085, normal

#4c: Model differences in dispersed seed diversity between mature & regenerating forests####
#Models summarized in Table S1 of article supplement

#lm used due to normal distribution of data (based on Shapiro tests)
#Models do not hold transect as random effect
summary(lm(richness~regeneration,dispersal_ren)) #p=0.613
summary(lm(shannons~regeneration,dispersal_ren)) #p=0.946
summary(lm(inv_simpson~regeneration,dispersal_ren)) #p=0.725
#Findings:
#Dispersed seed diversity does not significantly between intact and regenerating forests.

#4d: Visualize mean dispersed seeds community metrics between forest types####

#richness
dispersal_richness_plot<-
  ggplot(dispersal_ren,aes(x=regeneration, y=richness,fill=regeneration))+
  geom_boxplot(show.legend=FALSE)+
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  scale_x_discrete(label = c("Mature", "Regenerating"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle=15,vjust=0.6))+
  theme(plot.background = element_rect(fill='transparent',color='transparent'))+
  labs(y = ("Richness"), x = "")+
annotation_custom(grid::textGrob(x = unit(0.5,"npc"), y = unit(0.9,"npc"), label = "NS", gp= gpar(cex=1.5)))+
  annotation_custom(grid::textGrob(x = unit(0.05,"npc"), y = unit(0.93,"npc"), label = "d", gp= gpar(cex=2)))

#Hill-Shannon diversity
dispersal_shannon_plot<-
  ggplot(dispersal_ren,aes(x=regeneration, y=shannons,fill=regeneration))+
  geom_boxplot(show.legend=FALSE)+
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  scale_x_discrete(label = c("Mature", "Regenerating"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle=15,vjust=0.6))+
  theme(plot.background = element_rect(fill='transparent',color='transparent'))+
  labs(y = ("Shannon Index"), x = "")+
  annotation_custom(grid::textGrob(x = unit(0.5,"npc"), y = unit(0.9,"npc"), label = "NS", gp= gpar(cex=1.5)))+
  annotation_custom(grid::textGrob(x = unit(0.05,"npc"), y = unit(0.93,"npc"), label = "e", gp= gpar(cex=2)))

#inverse Simpson diversity
dispersal_evenness_plot<-
  ggplot(dispersal_ren,aes(x=regeneration, y=inv_simpson,fill=regeneration))+
  geom_boxplot(show.legend=FALSE)+
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  scale_x_discrete(label = c("Mature", "Regenerating"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle=15,vjust=0.6))+
  theme(plot.background = element_rect(fill='transparent',color='transparent'))+
  labs(y = ("Inverse Simpson Index"), x = "")+
  annotation_custom(grid::textGrob(x = unit(0.5,"npc"), y = unit(0.9,"npc"), label = "NS", gp= gpar(cex=1.5)))+
  annotation_custom(grid::textGrob(x = unit(0.05,"npc"), y = unit(0.93,"npc"), label = "f", gp= gpar(cex=2)))
  
#Figure S3 d,e,f in article supplement
dispersed.diversity<- grid.arrange(dispersal_richness_plot,dispersal_shannon_plot, dispersal_evenness_plot,
                                   ncol=3,
                               top=textGrob("Dispersed Seeds",gp=gpar(fontsize=20,font=2)))

#Figure S3 in article supplement
#Combine richness, Hill-Shannon, & inverse simpson from dispersed seeds with total seed rain figures (Part 3d) 
grid.arrange(rain.diversity,dispersed.diversity,ncol=1)


#4e: Model differences in dispersed seed abundance between forest types & over time####
#Models summarized in Table 2 in article and Table S2 in supplement

#testing for normality
shapiro.test(rain.dispersal_rates$dispersal.encounter.rate) #p=4.533e-9 rates non-normal
hist(rain.dispersal_rates$dispersal.encounter.rate)

shapiro.test(rain.dispersal_rates$nonnative.dispersal.encounter.rate) #p=4.655e-10, rates non-normal
hist(rain.dispersal_rates$nonnative.dispersal.encounter.rate)

#GLMMs chosen to account generalize for non-normal data and set transect as random effect.
#Tweedie distribution selected for positive continuous proportions with higher number of zeros

#Total dispersal models


#total dispersal rates without fruiting
summary(glmmTMB(dispersal.encounter.rate~regeneration+survey.number+(1|transect),
                family=tweedie,data=rain.dispersal_rates)) #forest type & week significant

#total dispersal rates with fruiting
summary(glmmTMB(dispersal.encounter.rate~regeneration+survey.number+fruiting.proportion+(1|transect),
                family=tweedie,data=rain.dispersal_rates)) #No sig explainitory variables

#Findings:
#Seed dispersal abundance decreases over time during our 7 week study.
#Regenerating forests have increased seed dispersal compared to mature forests
#Effects of forest type may be partially explained by higher proportions fo fruiting plants in regenerating forests, 
#but fruiting and forest type do not explain dispersal rates when considered together.


#rain.dispersal_rates$regeneration=as.factor(rain.dispersal_rates$regeneration)
#test.model<-glmmTMB(dispersal.encounter.rate~regeneration+survey.number+(1|transect),
              #  family=tweedie,data=rain.dispersal_rates) #forest type & week significant

#test.predict <- ggpredict(test.model,terms=c("survey.number","regeneration"))


#ggplot(test.predict,aes(as.factor(x),predicted))+geom_boxplot()+
 # geom_line(data=rain.dispersal_rates,x=as.factor(rain.dispersal_rates$survey.number),
   #         y=rain.dispersal_rates$dispersal.encounter.rate,color=rain.dispersal_rates$regeneration)
#plot(test.predict)



#Non-native dispersal models

#non-native dispersal rates without fruiting
summary(glmmTMB(nonnative.dispersal.encounter.rate~regeneration+survey.number+(1|transect),
                family=tweedie,data=rain.dispersal_rates)) ##forest type & week significant

#non-native dispersal rates with fruiting
summary(glmmTMB(nonnative.dispersal.encounter.rate~regeneration+survey.number+nonnative.fruiting.proportion+(1|transect),
                family=tweedie,data=rain.dispersal_rates))  #forest & type week significant, nonnative fruiting NOT sig

#Findings:
#Non-native seed dispersal abundance decreases over time during our 7 week study.
#Regenerating forests have increased non-native seed dispersal compared to mature forests
#Higher proportions of non-native fruiting plants in regenerating forests does not explain increased non-native seed dispersal  

non.native_dispersal_proportion <- (rain.dispersal_rates$nonnative.dispersal.count/rain.dispersal_rates$dispersal.count)
      
rain.dispersal_rates$non.native_dispersal_proportion <- non.native_dispersal_proportion

kruskal.test(non.native_dispersal_proportion~regeneration, data=rain.dispersal_rates) #p=0.0004829 sig difference in non-native proportion of seeds between habitats


#4f: Visualize dispersed seed abundance between forest types & over time####

total_dispersal_plot <- 
ggplot(rain.dispersal_rates,aes(x=as.factor(survey.number), dispersal.encounter.rate,fill=regeneration))+
  # geom_point()+
  geom_smooth(aes(group=regeneration,color=regeneration),method="glm",show.legend = FALSE)+
  geom_boxplot(show.legend=FALSE)+ 
  scale_fill_manual(values = c("slateblue", "skyblue"))+
  scale_color_manual(values = c("slateblue", "skyblue"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=16,face="bold"),
        plot.margin = unit(c(0,0,1.5,1.5), "cm"),
        axis.title.x=element_blank(),
  )+ ylim(0,3)+
  labs(y = ("Total Seed Rain\n(Seeds per Sampling Effort)"))+
  #geom_smooth(aes(group=Regeneration,color=Regeneration),method="lm",show.legend = FALSE)
  annotate(geom = "text", x = 0.75, y=17, label = "c", cex = 8)



non.native_dispersal_plot <- 
  ggplot(rain.dispersal_rates,aes(x=as.factor(survey.number), y=nonnative.dispersal.encounter.rate,fill=regeneration))+
  geom_smooth(aes(group=regeneration,color=regeneration),method="glm",show.legend = FALSE)+
  geom_boxplot(show.legend=FALSE)+
  #geom_boxplot(show.legend=FALSE)+ 
  scale_fill_manual(values = c("slateblue", "skyblue"))+
  scale_color_manual(values = c("slateblue", "skyblue"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=16,face="bold"),
        plot.margin = unit(c(0,0,1.5,1.5), "cm"),
        axis.title.x=element_blank(),
  )+ ylim(0,3)+
  labs(y = ("Non-native Dispersed Seeds\n(Seeds per Sampling Effort)"))+
  #geom_smooth(aes(group=Regeneration,color=Regeneration),method="lm",show.legend = FALSE)
  annotate(geom = "text", x = 0.75, y =2.9, label = "f", cex = 8)

#Figure 4 e,f in article
dispersed.plots<- grid.arrange(total_dispersal_plot,non.native_dispersal_plot,ncol=2,
             top=textGrob("Dispersed Seeds",gp=gpar(fontsize=20,font=2),x=0,hjust=0,vjust=-0.2))


#Figure 4 in article
#combine fruiting, seed rain, & dispersal plots
grid.arrange(fruit.plots,rain.plots,dispersed.plots,ncol=1,
             bottom=textGrob("Weeks",gp=gpar(fontsize=18,font=2)),
             top=textGrob("Fruiting Plants",gp=gpar(fontsize=20,font=2)))
#-------------------------------------------------------------------------------------------------------------------------------------------------------####
#PART 5: Seed Size ####

#5a: load & configure data ####
rain_individuals <- read.csv('Final Data/rain.seed.individuals.csv')
dispersal_individuals <- read.csv('Final Data/dispersal.seed.individuals.csv')

#5b: Compare differences in seed sizes between habitats ####

#check for normality
shapiro.test(rain_individuals$length.avg) #p<2.2e-16, non-normal
hist(rain_individuals$length.avg)

shapiro.test(rain_individuals$width.avg) #p<2.2e-16, non-normal
hist(rain_individuals$width.avg)
 
shapiro.test(dispersal_individuals$length.avg) #p=2.954e-13, non-normal
hist(dispersal_individuals$length.avg)

shapiro.test(dispersal_individuals$width.avg) #p=2.677e-15, non-normal
hist(dispersal_individuals$width.avg)

#Compare using non-parametric test
kruskal.test(length.avg~regeneration, data=rain_individuals) #p=2.09e-11
kruskal.test(width.avg~regeneration, data=rain_individuals) #p=2.2e-16
kruskal.test(length.avg~regeneration, data=dispersal_individuals) #p=1.26e-8
kruskal.test(width.avg~regeneration, data=dispersal_individuals) #p=2.73e-8
#Findings: Regenerating forests have shorter and thinner seeds in both rain and dispersal

#5c: Visualize seed sizes between forest types ####

#Plot seed rain length
rain.length.plot <- ggplot(rain_individuals,aes(x=regeneration, y=length.avg,fill=regeneration))+
  geom_boxplot(show.legend=FALSE)+
  #geom_violin()+
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  theme_bw()+
  ylim(0,25)+
  scale_x_discrete(label = c("Mature", "Regenerating"))+
theme(axis.text = element_text(size=14,color="black"),
      axis.title = element_text(size=20,face="bold"),
      axis.text.x = element_text(angle=15,vjust=0.6))+
  theme(plot.background = element_rect(fill='transparent',color='transparent'))+
  labs(y = ("Seed Length (mm)"), x = "")+
  annotation_custom(grid::textGrob(x = unit(0.5,"npc"), y = unit(0.9,"npc"), label = "***", gp= gpar(cex=1.5)))+
  annotation_custom(grid::textGrob(x = unit(0.05,"npc"), y = unit(0.93,"npc"), label = "a", gp= gpar(cex=2)))

#Plot seed rain width
rain.width.plot <- ggplot(rain_individuals,aes(x=regeneration, y=width.avg,fill=regeneration))+
  geom_boxplot(show.legend=FALSE)+
  #geom_violin()+
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  theme_bw()+
  ylim(0,25)+
  scale_x_discrete(label = c("Mature", "Regenerating"))+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle=15,vjust=0.6))+
  theme(plot.background = element_rect(fill='transparent',color='transparent'))+
  labs(y = ("Seed Width (mm)"), x = "")+
  annotation_custom(grid::textGrob(x = unit(0.5,"npc"), y = unit(0.9,"npc"), label = "***", gp= gpar(cex=1.5)))+
  annotation_custom(grid::textGrob(x = unit(0.05,"npc"), y = unit(0.93,"npc"), label = "b", gp= gpar(cex=2)))

rain.size.plots<- grid.arrange(rain.length.plot,rain.width.plot,ncol=2,
                               top=textGrob("Seed Rain",gp=gpar(fontsize=20,font=2),x=0,hjust=0,vjust=-0.2))

#Plot dispersed seed length
dispersed.length.plot <- ggplot(dispersal_individuals,aes(x=regeneration, y=length.avg,fill=regeneration))+
  geom_boxplot(show.legend=FALSE)+
  #geom_violin()+
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  theme_bw()+
  ylim(0,25)+
  scale_x_discrete(label = c("Mature", "Regenerating"))+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle=15,vjust=0.6))+
  theme(plot.background = element_rect(fill='transparent',color='transparent'))+
  labs(y = ("Seed Length (mm)"), x = "")+
  annotation_custom(grid::textGrob(x = unit(0.5,"npc"), y = unit(0.9,"npc"), label = "***", gp= gpar(cex=1.5)))+
  annotation_custom(grid::textGrob(x = unit(0.05,"npc"), y = unit(0.93,"npc"), label = "c", gp= gpar(cex=2)))


#Plot dispersed seed width
dispersed.width.plot <- ggplot(dispersal_individuals,aes(x=regeneration, y=width.avg,fill=regeneration))+
  geom_boxplot(show.legend=FALSE)+
  #geom_violin()+
  scale_fill_manual(values = (c("slateblue", "skyblue")))+
  theme_bw()+
  ylim(0,25)+
  scale_x_discrete(label = c("Mature", "Regenerating"))+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(angle=15,vjust=0.6))+
  theme(plot.background = element_rect(fill='transparent',color='transparent'))+
  labs(y = ("Seed Width (mm)"), x = "")+
  annotation_custom(grid::textGrob(x = unit(0.5,"npc"), y = unit(0.9,"npc"), label = "***", gp= gpar(cex=1.5)))+
  annotation_custom(grid::textGrob(x = unit(0.05,"npc"), y = unit(0.93,"npc"), label = "d", gp= gpar(cex=2)))

dispersed.size.plots<- grid.arrange(dispersed.length.plot,dispersed.width.plot,ncol=2,
                                    bottom=textGrob("Forest Type",gp=gpar(fontsize=18,font=2)),
                                    top=textGrob("Dispersed Seeds",gp=gpar(fontsize=20,font=2),x=0,hjust=0,vjust=-0.05))

grid.arrange(rain.size.plots,dispersed.size.plots,ncol=1,
             top=textGrob("",gp=gpar(fontsize=20,font=2)))








rain.dispersal_rates$regeneration=as.factor(rain.dispersal_rates$regeneration)
test.model<-glmmTMB(dispersal.encounter.rate~regeneration+survey.number+(1|transect),
                    family=tweedie,data=rain.dispersal_rates) #forest type & week significant

test.predict <- ggpredict(test.model,terms=c("survey.number","regeneration"))


ggplot(test.predict,aes(as.factor(x[]),predicted))+geom_boxplot()
geom_line(data=rain.dispersal_rates,x=as.factor(survey.number),y=dispersal.encounter.rate,color=regeneration)
plot(test.predict)


ggplot(data=test.predict,aes(x=survey.number,y=predicted))+geom_line()+geom_ribbon(data=test.predict,aes(ymin=conf.low,ymax=conf.high))
plot(test.predict)

total_dispersal_plot <- 
  ggplot(rain.dispersal_rates,aes(x=as.factor(survey.number), dispersal.encounter.rate,fill=regeneration))+
  # geom_point()+
  geom_ribbon(data=test.predict,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high,fill=group))+
  geom_smooth(aes(group=regeneration,color=regeneration),show.legend = FALSE)+ #,method="glm"
  geom_boxplot(show.legend=FALSE)+ 
  scale_fill_manual(values = c("slateblue", "skyblue"))+
  scale_color_manual(values = c("slateblue", "skyblue"))+
  theme_bw()+
  theme(axis.text = element_text(size=14,color="black"),
        axis.title = element_text(size=16,face="bold"),
        plot.margin = unit(c(0,0,1.5,1.5), "cm"),
        axis.title.x=element_blank(),
  )+ ylim(0,3)+
  labs(y = ("Total Dispersed Seed\n(Seeds per Sampling Effort)"))+
  #geom_smooth(aes(group=Regeneration,color=Regeneration),method="lm",show.legend = FALSE)
  annotate(geom = "text", x = 0.75, y =2.9, label = "e", cex = 8)
