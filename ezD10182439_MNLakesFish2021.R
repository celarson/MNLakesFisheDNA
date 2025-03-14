#Coding for MN Lakes Fish analysis

#Packages
library(ggplot2)
library(forcats)
library(reshape)
library(RColorBrewer)
library(tidyverse)
library(MetBrewer)
library(BiodiversityR) # also loads vegan
library(ggsci)
library(vegan)
library(fossil)
library(untb)

#Programs
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
  library(doBy)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # Collapse the data
  formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
  datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  
  # Rename columns
  names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
  names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
  names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


#Color vectors
eDNASampleTypesNoChoiceColVec<-met.brewer(name="Austria",n=7)
LakeClassColVec<-rev(met.brewer(name="Greek",n=3))
LakeColVec<-c("#A85E1A","#A19482","#7397BC","#045a8d","navy","#99d8c9","#238b45","#0E433A")
SmallLakeColVec<-c("#A85E1A","#A19482")
MedLakeColVec<-c("#7397BC","#045a8d","navy")
LargeLakeColVec<-c("#99d8c9","#238b45","#0E433A")
Habitatcolvec<-met.brewer(name="Juarez",n=3)
################
#Upload datasets
############

#Fish sequence abundance
MNFishSeq<-read.csv("MNDNRFishSequenceAbundance_UpdatedDec2024.csv",header=T)
#Change order of lake variable so smallest to largest
MNFishSeq$Lake<- factor(MNFishSeq$Lake, levels = c("East Chub", "Dunnigan", "Little Pine", "Loon", "Clark", "Big Pine", "Ripple", "Ely"))
#Reverse order of species so alphabetical
MNFishSeq$CommonName<-fct_rev(MNFishSeq$CommonName)
#Specify order of habitats
MNFishSeq$Habitat<- factor(MNFishSeq$Habitat, levels = c("deep point", "open water", "nearshore", "fringe", "inlet", "outlet", "ramp", "choice"))
#reverse order of depth
MNFishSeq$SampDepth<-fct_rev(MNFishSeq$SampDepth)

#Clean dataset
#how many reads to controls have
subset(MNFishSeq, Group=="No match")
#1389 reads no match 1318b3 Loon open water bottom AM12S
subset(MNFishSeq, Group=="NO RESULTS")
#7 samples 1315s,1358s, 1400s, 1411s, 1414s, 1539s, 1611m3
subset(MNFishSeq, Group=="control fish")
#None
subset(MNFishSeq, Group=="CA fish")
#None
subset(MNFishSeq,Group=="AIS")
#Ripple 1430-s1 tubenose goby (5 reads) and ripple 1430-s2 Round goby (2468 reads)
#subset those samples to see what other species show up
subset(MNFishSeq,SampID=="1430")
ggplot(subset(MNFishSeq,SampID=="1430"), aes(fill=CommonName, y=SignalReads, x=marker)) + 
  geom_bar(position="stack", stat="identity", color="black")+
  scale_fill_brewer(palette = "Paired")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  facet_wrap(.~SampDepthRepID)
  
#Exclude "NO RESUlTS"
MNFishSeqComplete<-subset(MNFishSeq, Group!="NO RESULTS")

#What is the range/mean/median signal reads
#12S
summary(subset(MNFishSeqComplete, marker=="Am12S")$SignalReads)
#Min=5    1stq=1978    median=6779   mean=11299   3rdq=16576  max=101683 
subset(MNFishSeqComplete, marker=="Am12S") %>%
  arrange(SignalReads) %>%  # arrange in descending order
  slice(1:20)                      # return rows 1 through 10

#16S
summary(subset(MNFishSeqComplete, marker=="Ac16S")$SignalReads)
#Min=10    1stq=959    median=3751   mean=6956   3rdq=9873  max=147648
subset(MNFishSeqComplete, marker=="Ac16S") %>%
  arrange(SignalReads) %>%  # arrange in descending order
  slice(1:20)                      # return rows 1 through 10



#delete non fish reads
MNFishSeqfilt<-subset(MNFishSeqComplete, Group!="nonfish")
#delete no match read
MNFishSeqNoContfish<-subset(MNFishSeqfilt, Group!="No match")


#Create site sequences rather than sample replicate
MNFishSeqSiteagg<- aggregate(MNFishSeqNoContfish$SignalReads, 
                             by=list(MNFishSeqNoContfish$Lake,
                                     MNFishSeqNoContfish$CommonName,
                                     MNFishSeqNoContfish$SampID,MNFishSeqNoContfish$Depth,
                                     MNFishSeqNoContfish$Habitat,
                                     MNFishSeqNoContfish$ZoneHabitat), FUN=sum)
#name variables
names(MNFishSeqSiteagg) <- c("Lake.Name", "Species.Common.Name", "SiteID","SampleDepth","Habitat","Zone","SeqAbundance")
#Add presence absence variable
MNFishSeqSiteagg$Presence<-ifelse(MNFishSeqSiteagg$SeqAbundance<1,0,1)

#Relative sequence variable
MNFishSeqRelSeq<-MNFishSeqSiteagg%>% 
  group_by(SiteID) %>% 
  mutate(relAbund = SeqAbundance / sum(SeqAbundance))

#Site information
MNFishSite<-read.csv("MNDNRSiteInfoMNLakes.csv")

#merge sequence aggregated dataset with site info
MNFishagg<-merge(MNFishSeqSiteagg, MNFishSite, by= c("Lake.Name","SiteID", "Habitat","SampleDepth"), all.x=T)
MNFishagg<-MNFishagg%>%mutate(Zone= cut(SiteDepth_m, breaks = c(-Inf,1.5,4.5,Inf)))


#DNR Dataset
DNRFish<-read.csv("ezD10182435_MNDNRTradData.csv",header=T)
#Add presence absence variable
DNRFish$Presence<-ifelse(DNRFish$Station.Total.Fish.Count<1,0,1)
#Remove Rock and Ripple
DNRFish<-subset(DNRFish,Lake.Name!="Ripple" & Lake.Name!="Rock")
#Change order of levels
DNRFish$Lake.Name<- factor(DNRFish$Lake.Name, levels = c("East Chub", "Dunnigan", "Little Pine", "Loon", "Clark", "Big Pine", "Ely"))
#Reverse order of species so alphabetical
DNRFish$Species.Common.Name<-fct_rev(DNRFish$Species.Common.Name)
#Only include the most recent survey
DNRFish$Surveys<-paste(DNRFish$Lake.Name,DNRFish$Survey.ID.Date,DNRFish$Station.Type.Abbreviation,DNRFish$Survey.Type.Abbreviation)

##############
#Coding
###############

#aggregate to number of samples required to get signal per lake
aggsumMNFishSeq<- aggregate(MNFishSeqRelSeq$Presence, by=list(MNFishSeqRelSeq$Lake.Name,MNFishSeqRelSeq$Species.Common.Name), FUN=sum)

#visualize with heatmap
ggplot(aggsumMNFishSeq, aes(Group.1, Group.2))+
  geom_tile(aes(fill = x))+
  scale_fill_gradient(low = "#ffffd9", high = "#081d58", name="Sum\n# Positive\nSamples")+
  labs(x="Lake", y="Common Name")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))

#Create relative
MNFishPresMatrix<- cast(MNFishSeqRelSeq, Lake.Name + SiteID + SampleDepth + Habitat + Zone ~ Species.Common.Name, value='Presence')
#replace NA's with 0's
MNFishPresMatrix[is.na(MNFishPresMatrix)] <- 0
#melt back so there are 0's
MNFishPres0s<-melt(MNFishPresMatrix, na.rm = FALSE, value.name ="Presence", id =c("Lake.Name","SiteID","Depth","Habitat","Zone"))

#aggregate to average number of samples required to get signal per lake
aggmeanMNFishPres0s<- aggregate(MNFishPres0s$value, by=list(MNFishPres0s$Lake.Name,MNFishPres0s$Species.Common.Name), FUN=mean)
#replace 0 with NA for visualization
aggmeanMNFishPres<-subset(aggmeanMNFishPres0s,x>0)
#visualize with heatmap
ggplot(aggmeanMNFishPres, aes(Group.1, Group.2))+
  geom_tile(aes(fill = x))+
  scale_fill_gradient(low = "#ffffd9", high = "#081d58", name="Mean\n# Positive\nSamples")+
  labs(x="Lake", y="Common Name")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))

#Species accumulation curve
#Create species by site matrix abundance
MNFishRelMatrix<-cast(MNFishSeqRelSeq, Lake.Name + SiteID + SampleDepth + Habitat + Zone ~ Species.Common.Name, value='relAbund')
MNFishRelMatrix<- as.data.frame(MNFishRelMatrix)
MNFishRelMatrix[is.na(MNFishRelMatrix)] <- 0

#split into community and environmental
MNFishRelMatrixComm<-MNFishRelMatrix[,6:ncol(MNFishRelMatrix)]
MNFishRelMatrixEnv<-MNFishRelMatrix[,1:5]

#Now do again for each lake
#Ely split into community and environmental
MNFishRelMatrixCommEly<-subset(MNFishRelMatrix, Lake.Name=="Ely")[,6:ncol(MNFishRelMatrix)]
MNFishRelMatrixCommEly<-MNFishRelMatrixCommEly[, which(colSums(MNFishRelMatrixCommEly) != 0)]
#26 species in Ely, 73 sites
MNFishRelMatrixEnvEly<-subset(MNFishRelMatrix, Lake.Name=="Ely")[,1:5]

#Ripple split into community and environmental
MNFishRelMatrixCommRip<-subset(MNFishRelMatrix, Lake.Name=="Ripple")[,6:ncol(MNFishRelMatrix)]
MNFishRelMatrixCommRip<-MNFishRelMatrixCommRip[, which(colSums(MNFishRelMatrixCommRip) != 0)]
#31 species in Ripple, 69 sites
MNFishRelMatrixEnvRip<-subset(MNFishRelMatrix, Lake.Name=="Ripple")[,1:5]

#Big Pine split into community and environmental
MNFishRelMatrixCommBP<-subset(MNFishRelMatrix, Lake.Name=="Big Pine")[,6:ncol(MNFishRelMatrix)]
MNFishRelMatrixCommBP<-MNFishRelMatrixCommBP[, which(colSums(MNFishRelMatrixCommBP) != 0)]
#25 species in Big Pine, 71 sites
MNFishRelMatrixEnvBP<-subset(MNFishRelMatrix, Lake.Name=="Big Pine")[,1:5]
MNFishRelMatrixEnvBP$Zone<-as.factor(MNFishRelMatrixEnvBP$Zone)

#no choice big pine
MNFishRelMatrixCommBPnc<-subset(MNFishRelMatrix, Lake.Name=="Big Pine" & Habitat!="choice")[,6:ncol(MNFishRelMatrix)]
MNFishRelMatrixCommBPnc<-MNFishRelMatrixCommBPnc[, which(colSums(MNFishRelMatrixCommBPnc) != 0)]
#25 species in Big Pine, 66 sites
MNFishRelMatrixEnvBPnc<-subset(MNFishRelMatrix, Lake.Name=="Big Pine"& Habitat!="choice")[,1:5]

#Clark split into community and environmental
MNFishRelMatrixCommCla<-subset(MNFishRelMatrix, Lake.Name=="Clark")[,6:ncol(MNFishRelMatrix)]
MNFishRelMatrixCommCla<-MNFishRelMatrixCommCla[, which(colSums(MNFishRelMatrixCommCla) != 0)]
#25 species in Clark, 45 sites
MNFishRelMatrixEnvCla<-subset(MNFishRelMatrix, Lake.Name=="Clark")[,1:5]

#Loon split into community and environmental
MNFishRelMatrixCommLoon<-subset(MNFishRelMatrix, Lake.Name=="Loon")[,6:ncol(MNFishRelMatrix)]
MNFishRelMatrixCommLoon<-MNFishRelMatrixCommLoon[, which(colSums(MNFishRelMatrixCommLoon) != 0)]
#23 species in Loon, 32 sites
MNFishRelMatrixEnvLoon<-subset(MNFishRelMatrix, Lake.Name=="Loon")[,1:5]

#Little Pine split into community and environmental
MNFishRelMatrixCommLP<-subset(MNFishRelMatrix, Lake.Name=="Little Pine")[,6:ncol(MNFishRelMatrix)]
MNFishRelMatrixCommLP<-MNFishRelMatrixCommLP[, which(colSums(MNFishRelMatrixCommLP) != 0)]
#28 species in Little Pine, 35 sites
MNFishRelMatrixEnvLP<-subset(MNFishRelMatrix, Lake.Name=="Little Pine")[,1:5]

#Dunnigan split into community and environmental
MNFishRelMatrixCommDun<-subset(MNFishRelMatrix, Lake.Name=="Dunnigan")[,6:ncol(MNFishRelMatrix)]
MNFishRelMatrixCommDun<-MNFishRelMatrixCommDun[, which(colSums(MNFishRelMatrixCommDun) != 0)]
#11 species in Dunnigan, 17 sites
MNFishRelMatrixEnvDun<-subset(MNFishRelMatrix, Lake.Name=="Dunnigan")[,1:5]

#East Chub split into community and environmental
MNFishRelMatrixCommEC<-subset(MNFishRelMatrix, Lake.Name=="East Chub")[,6:ncol(MNFishRelMatrix)]
MNFishRelMatrixCommEC<-MNFishRelMatrixCommEC[, which(colSums(MNFishRelMatrixCommEC) != 0)]
#9 species in East Chub, 16 sites
MNFishRelMatrixEnvEC<-subset(MNFishRelMatrix, Lake.Name=="East Chub")[,1:5]
MNFishRelMatrixEnvEC$Zone<-as.factor(MNFishRelMatrixEnvEC$Zone)

#species accumulation curve plotting
MNFishPoolAccum<-poolaccum(MNFishRelMatrixComm)
plot(MNFishPoolAccum)
MNFishaccum<-specaccum(MNFishRelMatrixComm, method = "exact", conditioned =TRUE, 
                         gamma = "jack1",  w = NULL)
plot(MNFishaccum)

MNFishAccum.1<-accumcomp(MNFishRelMatrixComm, y=MNFishRelMatrixEnv, factor='Lake.Name', 
                     method='exact', conditioned=FALSE, plotit=FALSE)
MNFishaccum.long1 <- accumcomp.long(MNFishAccum.1, ci=NA, label.freq=5)
MNFishaccum.long1$Grouping<- factor(MNFishaccum.long1$Grouping, levels = c("East Chub", "Dunnigan", "Little Pine", "Loon", "Clark", "Big Pine", "Ripple", "Ely"))

ggplot(MNFishaccum.long1, aes(x = Sites, y = Richness, colour=Grouping)) + 
  geom_line(linewidth=2) +
  labs(x = "Number of samples", y = "Number of Species",
       colour="Lake Name")+
  scale_color_manual(values=c(LakeColVec))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),legend.title=element_text(size=14))

ggplot(subset(MNFishaccum.long1, Grouping=="Ely"), aes(x = Sites, y = Richness)) + 
  geom_line(linewidth=2) +
  labs(x = "Number of samples", y = "Number of Fish Species")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),legend.title=element_text(size=14))


#How many samples does it take to get 95% of species richness, without taking into account habitat
#Total richness for each lake
chao1(MNFishRelMatrixCommEly,taxa.row=F)
#Ely ch1 26, 95% is 24.7, and 50% is 13
chao2(MNFishRelMatrixCommEly,taxa.row=F)
#Ely ch2 28, 95% is 26.6, and 50% is 14

chao1(MNFishRelMatrixCommRip,taxa.row=F)
#Ripple ch1 31, 95% is 29.45, and 50% is 15.5
chao2(MNFishRelMatrixCommRip,taxa.row=F)
#Ripple ch2 43.5, 95% is 41.325, and 50% is 21.75

chao1(MNFishRelMatrixCommBP,taxa.row=F)
#Big ch1 Pine 25, 95% is 23.75, and 50% is 12.5
chao2(MNFishRelMatrixCommBP,taxa.row=F)
#Big ch2 Pine 25, 95% is 23.75, and 50% is 12.5

chao1(MNFishRelMatrixCommCla,taxa.row=F)
#Clark ch1 25, 95% is 23.75, and 50% is 12.5
chao2(MNFishRelMatrixCommCla,taxa.row=F)
#Clark ch2 35, 95% is 23.75, and 50% is 12.5

chao1(MNFishRelMatrixCommLoon,taxa.row=F)
#Loon ch1 23, 95% is 21.85, and 50% is 11.5
chao2(MNFishRelMatrixCommLoon,taxa.row=F)
#Loon ch2 34, 95% is 22.8, and 50% is 17

chao1(MNFishRelMatrixCommLP,taxa.row=F)
#Little Pine ch1 28, 95% is 26.6, and 50% is 14
chao2(MNFishRelMatrixCommLP,taxa.row=F)
#Little Pine ch2 32.5, 95% is 30.874, and 50% is 16.25

chao1(MNFishRelMatrixCommDun,taxa.row=F)
#Dunnigan ch1 11, 95% is 10.45, and 50% is 5.5
chao2(MNFishRelMatrixCommDun,taxa.row=F)
#Dunnigan ch2 13.67, 95% is 12.99, and 50% is 6.84

chao1(MNFishRelMatrixCommEC,taxa.row=F)
#East Chub ch1 9, 95% is 8.55, and 50% is 4.5
chao2(MNFishRelMatrixCommEC,taxa.row=F)
#EAst Chub ch2 9, 95% is 8.55, and 50% is 4.5

#Plot each lake group separately
ggplot(subset(MNFishaccum.long1, Grouping=="Ely" | 
                Grouping=="Ripple" | Grouping=="Big Pine"), 
       aes(x = Sites, y = Richness, colour=Grouping)) + 
  geom_line(size=2) +
  scale_color_manual(values=c(LargeLakeColVec))+
  labs(x = "Number of samples", y = "Number of Fish Species",
       colour="Lake Name")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),legend.title=element_text(size=14))

ggplot(subset(MNFishaccum.long1, Grouping=="Little Pine" | 
                Grouping=="Loon" | Grouping=="Clark"), 
       aes(x = Sites, y = Richness, colour=Grouping)) + 
  geom_line(size=2) +
  scale_color_manual(values=c(MedLakeColVec))+
  labs(x = "Number of samples", y = "Number of Fish Species",
       colour="Lake Name")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),legend.title=element_text(size=14))

ggplot(subset(MNFishaccum.long1, Grouping=="East Chub" | 
                Grouping=="Dunnigan"), 
       aes(x = Sites, y = Richness, colour=Grouping)) + 
  geom_line(size=2) +
  scale_color_manual(values=c(SmallLakeColVec))+
  labs(x = "Number of samples", y = "Number of Fish Species",
       colour="Lake Name")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),legend.title=element_text(size=14))

#East chub and big pine are the most well sampled, so use those two for further analysis

MNFishAccum.BPH<-accumcomp(MNFishRelMatrixCommBP, y=MNFishRelMatrixEnvBP, factor='Zone', 
                         method='exact', conditioned=FALSE, plotit=FALSE)
MNFishaccum.longBPH<- accumcomp.long(MNFishAccum.BPH, ci=NA, label.freq=5)

ggplot(MNFishaccum.longBPH, aes(x = Sites, y = Richness, colour=Grouping)) + 
  geom_line(size=2) +
  labs(x = "Number of samples", y = "Number of Species",
       colour="Habitat")+
  scale_color_manual(values=c(Habitatcolvec))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),legend.title=element_text(size=14))

#NMDS Big Pine
BP_NMDS<-metaMDS(MNFishRelMatrixCommBP, distance = "bray", k = 3, 
                      maxit=1000, trymax = 300, wascores = FALSE,
                      autotransform = FALSE, trace = 2, noshare = FALSE)
#stress = 0.13

#Stressplot 
stressplot(BP_NMDS)

#NMDS plot for zone
ordiplot(BP_NMDS, type="n")
with(BP_NMDS, points(BP_NMDS, display="sites", col=Habitatcolvec[MNFishRelMatrixEnvBP$Zone], pch=19))
with(BP_NMDS, legend("topleft", legend=levels(MNFishRelMatrixEnvBP$Zone), bty="n", col=Habitatcolvec, pch=19, pt.bg=Habitatcolvec))
with(BP_NMDS, ordiellipse(BP_NMDS, MNFishRelMatrixEnvBP$Zone, kind="se", conf=0.95, lwd=2, col="#a82203", show.groups = "fringe"))
with(BP_NMDS, ordiellipse(BP_NMDS, MNFishRelMatrixEnvBP$Zone, kind="se", conf=0.95, lwd=2, col="#208cc0", show.groups = "littoral"))
with(BP_NMDS, ordiellipse(BP_NMDS, MNFishRelMatrixEnvBP$Zone, kind="se", conf=0.95, lwd=2, col="#f1af3a", show.groups = "open water"))

#Fringe vs. open water and littoral
#Triangle plot fringe vs. off shore

#Big Pine
#Subset so comm and env together
MNFishRelMatrixBP<-subset(MNFishRelMatrix, Lake.Name=="Big Pine")

#SUBSET SITES ACCORDING TO GEAR TYPE
fringeBP<-subset(MNFishRelMatrixBP, Zone=="fringe")
#28 obs
littoralBP<-subset(MNFishRelMatrixBP, Zone=="littoral")
#17 obs
openwaterBP<-subset(MNFishRelMatrixBP, Zone=="open water")
#26 obs
#PREPARE BLANK DATA FRAME FOR RANDOM DRAWS USING SAME METHODS AS TREBITZ ET AL. (2009) AND HOFFMAN ET AL. (2016)
n.samples<-expand.grid(n.fringe=seq(0,10,1), n.littoral=seq(0,10,1), n.openwater=seq(0,10,1))
n.samples$sum<-rowSums(n.samples)
n.samples10BP<-subset(n.samples, sum==10)
n.samples10BP$sum<-NULL
n.samples10BP$result<-0
BP100<-c(1:100)
#RUN LOOPS
for (i in 1:nrow(n.samples10BP))
{
  for (j in 1:length(BP100))
  {
    fringe.sample<-fringeBP[sample(1:nrow(fringeBP), n.samples10BP[i,1], replace=FALSE),]
    littoral.sample<-littoralBP[sample(1:nrow(littoralBP), n.samples10BP[i,2], replace=FALSE),]
    openwater.sample<-openwaterBP[sample(1:nrow(openwaterBP), n.samples10BP[i,3], replace=FALSE),]
    sample<-rbind(fringe.sample, littoral.sample, openwater.sample)
    totals<-as.numeric(colSums(sample[,6:44]))
    BP100[j]<-specnumber(totals)
  }
  n.samples10BP[i,4]<-mean(BP100)
}

#find maximum
max(n.samples10BP$result)
#21.53
#GRAPH OF RANDOMIZATION RESULTS
ggplot(n.samples10BP, aes(n.fringe, n.littoral))+
  geom_tile(aes(fill = result), color='white')+
  scale_fill_gradient(low="lightgreen",high="black")+
  labs(x="Number of Fringe Samples",y="Number of Littoral Samples",
       fill="Fish Species\nRichness")+
  scale_y_continuous(breaks=c(0,2,4,6,8,10))+
  scale_x_continuous(breaks=c(0,2,4,6,8,10))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),legend.title=element_text(size=14))

#what is the rarest species in big pine
sort(colSums(MNFishRelMatrixBP[,6:44]))
#Blacknose Shiner, Blackchin Shiner, Bowfin, Bluntnose Minnow

#Try again with more samples
#PREPARE BLANK DATA FRAME FOR RANDOM DRAWS USING SAME METHODS AS TREBITZ ET AL. (2009) AND HOFFMAN ET AL. (2016)
n.samples20<-expand.grid(n.fringe=seq(0,20,1), n.littoral=seq(0,15,1), n.openwater=seq(0,20,1))
n.samples20$sum<-rowSums(n.samples20)
n.samples20BP<-subset(n.samples20, sum==20)
n.samples20BP$sum<-NULL
n.samples20BP$result<-0
BP100<-c(1:100)
#RUN LOOPS
for (i in 1:nrow(n.samples20BP))
{
  for (j in 1:length(BP100))
  {
    fringe.sample<-fringeBP[sample(1:nrow(fringeBP), n.samples20BP[i,1], replace=FALSE),]
    littoral.sample<-littoralBP[sample(1:nrow(littoralBP), n.samples20BP[i,2], replace=FALSE),]
    openwater.sample<-openwaterBP[sample(1:nrow(openwaterBP), n.samples20BP[i,3], replace=FALSE),]
    sample<-rbind(fringe.sample, littoral.sample, openwater.sample)
    totals<-as.numeric(colSums(sample[,6:44]))
    BP100[j]<-specnumber(totals)
  }
  n.samples20BP[i,4]<-mean(BP100)
}

#find maximum
max(n.samples20BP$result)
#24.02
#GRAPH OF RANDOMIZATION RESULTS
ggplot(n.samples20BP, aes(n.fringe, n.openwater))+
  geom_tile(aes(fill = result), color='white')+
  scale_fill_gradient(low="lightgreen",high="black")+
  labs(x="Number of Fringe Samples",y="Number of Open Water Samples",
       fill="Fish Species\nRichness")+
  scale_y_continuous(breaks=c(0,4,8,12,16,20))+
  scale_x_continuous(breaks=c(0,4,8,12,16,20))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),legend.title=element_text(size=14))


#EastChub
MNFishAccum.ECH<-accumcomp(MNFishRelMatrixCommEC, y=MNFishRelMatrixEnvEC, factor='Zone', 
                           method='exact', conditioned=FALSE, plotit=FALSE)
MNFishaccum.longECH<- accumcomp.long(MNFishAccum.ECH, ci=NA, label.freq=5)

ggplot(MNFishaccum.longECH, aes(x = Sites, y = Richness, colour=Grouping)) + 
  geom_line(size=2) +
  labs(x = "Number of samples", y = "Number of Species",
       colour="Habitat")+
  scale_color_manual(values=c(Habitatcolvec))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),legend.title=element_text(size=14))

#NMDS East Chub
EC_NMDS<-metaMDS(MNFishRelMatrixCommEC, distance = "bray", k = 3, 
                 maxit=1000, trymax = 300, wascores = FALSE,
                 autotransform = FALSE, trace = 2, noshare = FALSE)
#stress = 0.05

#Stressplot 
stressplot(EC_NMDS)

#NMDS plot for zone
ordiplot(EC_NMDS, type="n")
with(EC_NMDS, points(BP_NMDS, display="sites", col=Habitatcolvec[MNFishRelMatrixEnvEC$Zone], pch=19))
with(EC_NMDS, legend("topleft", legend=levels(MNFishRelMatrixEnvEC$Zone), bty="n", col=Habitatcolvec, pch=19, pt.bg=Habitatcolvec))
with(EC_NMDS, ordiellipse(EC_NMDS, MNFishRelMatrixEnvEC$Zone, kind="se", conf=0.95, lwd=2, col="#a82203", show.groups = "fringe"))
with(EC_NMDS, ordiellipse(EC_NMDS, MNFishRelMatrixEnvEC$Zone, kind="se", conf=0.95, lwd=2, col="#208cc0", show.groups = "littoral"))
with(EC_NMDS, ordiellipse(EC_NMDS, MNFishRelMatrixEnvEC$Zone, kind="se", conf=0.95, lwd=2, col="#f1af3a", show.groups = "open water"))

#Fringe vs. open water and littoral
#Triangle plot fringe vs. off shore

#East Chub
#Subset so comm and env together
MNFishRelMatrixEC<-subset(MNFishRelMatrix, Lake.Name=="East Chub")

#SUBSET SITES ACCORDING TO GEAR TYPE
fringeEC<-subset(MNFishRelMatrixEC, Zone=="fringe")
#6 obs
littoralEC<-subset(MNFishRelMatrixEC, Zone=="littoral")
#6 obs
openwaterEC<-subset(MNFishRelMatrixEC, Zone=="open water")
#4 obs
#PREPARE BLANK DATA FRAME FOR RANDOM DRAWS USING SAME METHODS AS TREBITZ ET AL. (2009) AND HOFFMAN ET AL. (2016)
n.samples<-expand.grid(n.fringe=seq(0,5,1), n.littoral=seq(0,5,1), n.openwater=seq(0,4,1))
n.samples$sum<-rowSums(n.samples)
n.samples5EC<-subset(n.samples, sum==5)
n.samples5EC$sum<-NULL
n.samples5EC$result<-0
EC100<-c(1:100)
#RUN LOOPS
for (i in 1:nrow(n.samples5EC))
{
  for (j in 1:length(EC100))
  {
    fringe.sample<-fringeEC[sample(1:nrow(fringeEC), n.samples5EC[i,1], replace=FALSE),]
    littoral.sample<-littoralEC[sample(1:nrow(littoralEC), n.samples5EC[i,2], replace=FALSE),]
    openwater.sample<-openwaterEC[sample(1:nrow(openwaterEC), n.samples5EC[i,3], replace=FALSE),]
    sample<-rbind(fringe.sample, littoral.sample, openwater.sample)
    totals<-as.numeric(colSums(sample[,6:44]))
    EC100[j]<-specnumber(totals)
  }
  n.samples5EC[i,4]<-mean(EC100)
}

#find maximum
max(n.samples5EC$result)
#8.92
#GRAPH OF RANDOMIZATION RESULTS
ggplot(n.samples5EC, aes(n.fringe, n.littoral))+
  geom_tile(aes(fill = result), color='white')+
  scale_fill_gradient(low="lightgreen",high="black")+
  labs(x="Number of Fringe Samples",y="Number of Littoral Samples",
       fill="Fish Species\nRichness")+
  scale_y_continuous(breaks=c(0,1,2,3,4,5))+
  scale_x_continuous(breaks=c(0,1,2,3,4,5))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),legend.title=element_text(size=14))

#####################
#Now see how different targetted approaches compare
###########

#Aggregate dataset based on sample type and lake
MNFishRelMatrixtargetagg<-aggregate(MNFishRelMatrix[,6:44], 
                                     by=list(MNFishRelMatrix$Lake,
                                             MNFishRelMatrix$Habitat), FUN=mean)
MNFishRelMatrixtargetagg$Richness<-species.count(MNFishRelMatrixtargetagg[,3:41])
#find each lake's richness
MNFishRelMatrixlakeagg<-aggregate(MNFishRelMatrix[,6:44], 
                                    by=list(MNFishRelMatrix$Lake), FUN=mean)
MNFishRelMatrixlakeagg$Richness<-species.count(MNFishRelMatrixlakeagg[,2:40])

MNFishRelMatrixtargetagg$TotalRichness<-MNFishRelMatrixtargetagg$Group.1
levels(MNFishRelMatrixtargetagg$TotalRichness)[levels(MNFishRelMatrixtargetagg$TotalRichness)=="East Chub"] <- 10
levels(MNFishRelMatrixtargetagg$TotalRichness)[levels(MNFishRelMatrixtargetagg$TotalRichness)=="Dunnigan"] <- 12
levels(MNFishRelMatrixtargetagg$TotalRichness)[levels(MNFishRelMatrixtargetagg$TotalRichness)=="Little Pine"] <- 29
levels(MNFishRelMatrixtargetagg$TotalRichness)[levels(MNFishRelMatrixtargetagg$TotalRichness)=="Loon"] <- 24
levels(MNFishRelMatrixtargetagg$TotalRichness)[levels(MNFishRelMatrixtargetagg$TotalRichness)=="Clark"] <- 26
levels(MNFishRelMatrixtargetagg$TotalRichness)[levels(MNFishRelMatrixtargetagg$TotalRichness)=="Big Pine"] <- 26
levels(MNFishRelMatrixtargetagg$TotalRichness)[levels(MNFishRelMatrixtargetagg$TotalRichness)=="Ripple"] <- 32
levels(MNFishRelMatrixtargetagg$TotalRichness)[levels(MNFishRelMatrixtargetagg$TotalRichness)=="Ely"] <- 27
MNFishRelMatrixtargetagg$TotalRichness<-as.numeric(as.character(MNFishRelMatrixtargetagg$TotalRichness))

MNFishRelMatrixtargetagg$PercSpecies<-MNFishRelMatrixtargetagg$Richness/MNFishRelMatrixtargetagg$TotalRichness
MNFishRelMatrixtargetagg$LakeSize<-MNFishRelMatrixtargetagg$Group.1
levels(MNFishRelMatrixtargetagg$LakeSize)[levels(MNFishRelMatrixtargetagg$LakeSize)=="East Chub"] <- "Small"
levels(MNFishRelMatrixtargetagg$LakeSize)[levels(MNFishRelMatrixtargetagg$LakeSize)=="Dunnigan"] <- "Small"
levels(MNFishRelMatrixtargetagg$LakeSize)[levels(MNFishRelMatrixtargetagg$LakeSize)=="Little Pine"] <- "Medium"
levels(MNFishRelMatrixtargetagg$LakeSize)[levels(MNFishRelMatrixtargetagg$LakeSize)=="Loon"] <- "Medium"
levels(MNFishRelMatrixtargetagg$LakeSize)[levels(MNFishRelMatrixtargetagg$LakeSize)=="Clark"] <- "Medium"
levels(MNFishRelMatrixtargetagg$LakeSize)[levels(MNFishRelMatrixtargetagg$LakeSize)=="Big Pine"] <- "Large"
levels(MNFishRelMatrixtargetagg$LakeSize)[levels(MNFishRelMatrixtargetagg$LakeSize)=="Ripple"] <- "Large"
levels(MNFishRelMatrixtargetagg$LakeSize)[levels(MNFishRelMatrixtargetagg$LakeSize)=="Ely"] <- "Large"

sumMNFishRelMatrixtargetagg<-summarySE(MNFishRelMatrixtargetagg, measurevar="PercSpecies",  groupvars=c("Group.2","LakeSize"))
sumMNFishRelMatrixtargetaggpic<-subset(sumMNFishRelMatrixtargetagg, Group.2=="ramp" | Group.2=="deep point" |
         Group.2=="inlet" | Group.2=="outlet")
ggplot(sumMNFishRelMatrixtargetaggpic, 
              aes(x=LakeSize, y=PercSpecies)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=PercSpecies-se, ymax=PercSpecies+se)) +  
  xlab("Lake Size")+
  ylab("Mean Percent Total Fish eDNA Richness")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),legend.title=element_text(size=14))+
  facet_wrap(.~Group.2)





#what is the rarest species in big pine
sort(colSums(MNFishRelMatrixECnc[,5:43]))
#Blacknose Shiner, Blackchin Shiner, Bowfin, Bluntnose Minnow
sort(colSums(MNFishRelMatrixBP[,5:43]))

#Try again with more samples
#PREPARE BLANK DATA FRAME FOR RANDOM DRAWS USING SAME METHODS AS TREBITZ ET AL. (2009) AND HOFFMAN ET AL. (2016)
n.samples20<-expand.grid(n.fringe=seq(0,20,1), n.littoral=seq(0,15,1), n.openwater=seq(0,20,1))
n.samples20$sum<-rowSums(n.samples20)
n.samples20BP<-subset(n.samples20, sum==20)
n.samples20BP$sum<-NULL
n.samples20BP$result<-0
BP100<-c(1:100)
#RUN LOOPS
for (i in 1:nrow(n.samples20BP))
{
  for (j in 1:length(BP100))
  {
    fringe.sample<-fringeBP[sample(1:nrow(fringeBP), n.samples20BP[i,1], replace=FALSE),]
    littoral.sample<-littoralBP[sample(1:nrow(littoralBP), n.samples20BP[i,2], replace=FALSE),]
    openwater.sample<-openwaterBP[sample(1:nrow(openwaterBP), n.samples20BP[i,3], replace=FALSE),]
    sample<-rbind(fringe.sample, littoral.sample, openwater.sample)
    totals<-as.numeric(colSums(sample[,5:43]))
    BP100[j]<-specnumber(totals)
  }
  n.samples20BP[i,4]<-mean(BP100)
}

#find maximum
max(n.samples20BP$result)
#24.1
#GRAPH OF RANDOMIZATION RESULTS
ggplot(n.samples20BP, aes(n.fringe, n.openwater))+
  geom_tile(aes(fill = result), color='white')+
  scale_fill_gradient(low="lightgreen",high="black")+
  labs(x="Number of Fringe samples",y="Number of Open Water samples",
       fill="Fish Species\nRichness")+
  scale_y_continuous(breaks=c(0,4,8,12,16,20))+
  scale_x_continuous(breaks=c(0,4,8,12,16,20))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),legend.title=element_text(size=14))






MNFishSeqRelSeqBPnc<-subset(MNFishSeqRelSeq, Lake.Name=="Big Pine" & Habitat!="choice")
MNFishSeqRelSeqBPnc$HabitatCondensed<-MNFishSeqRelSeqBPnc$Habitat
levels(MNFishSeqRelSeqBPnc$HabitatCondensed)[levels(MNFishSeqRelSeqBPnc$HabitatCondensed)=="deep point"] <- "open water"
levels(MNFishSeqRelSeqBPnc$HabitatCondensed)[levels(MNFishSeqRelSeqBPnc$HabitatCondensed)=="inlet"] <- "fringe"
levels(MNFishSeqRelSeqBPnc$HabitatCondensed)[levels(MNFishSeqRelSeqBPnc$HabitatCondensed)=="outlet"] <- "fringe"
levels(MNFishSeqRelSeqBPnc$HabitatCondensed)[levels(MNFishSeqRelSeqBPnc$HabitatCondensed)=="ramp"] <- "fringe"
levels(MNFishSeqRelSeqBPnc$HabitatCondensed)[levels(MNFishSeqRelSeqBPnc$HabitatCondensed)=="nearshore"] <- "littoral"
MNFishSeqRelSeqBPnc$HabitatCondensed<-droplevels(MNFishSeqRelSeqBPnc$HabitatCondensed)
#Condense into fringe vs. littoral and open water
MNFishSeqRelSeqBPnc$HabitatCondensed2<-MNFishSeqRelSeqBPnc$HabitatCondensed
levels(MNFishSeqRelSeqBPnc$HabitatCondensed2)[levels(MNFishSeqRelSeqBPnc$HabitatCondensed2)=="open water"] <- "off shore"
levels(MNFishSeqRelSeqBPnc$HabitatCondensed2)[levels(MNFishSeqRelSeqBPnc$HabitatCondensed2)=="littoral"] <- "off shore"
MNFishSeqRelSeqBPnc$HabitatCondensed2<-droplevels(MNFishSeqRelSeqBPnc$HabitatCondensed2)
MNFishSeqRelSeqBPncpres<-#####delete extra columns from dataset######
#SUBSET SITES ACCORDING TO GEAR TYPE
fringe<-subset(MNFishSeqRelSeqBPnc, HabitatCondensed2=="fringe")
#28 obs
offshore<-subset(MNFishSeqRelSeqBPnc, HabitatCondensed2=="off shore")
#38
#PREPARE BLANK DATA FRAME FOR RANDOM DRAWS USING SAME METHODS AS TREBITZ ET AL. (2009) AND HOFFMAN ET AL. (2016)
df11<-expand.grid(n.fringe=seq(0,10,1), n.offshore=seq(0,10,1))
df11$sum<-rowSums(df11)
df12<-subset(df11, sum==10)
df12$sum<-NULL
df12$result<-0
vec1<-c(1:100)
#RUN LOOPS
for (i in 1:nrow(df12))
{
  for (j in 1:length(vec1))
  {
    fringe.sample<-fringe[sample(1:nrow(fringe), df12[i,1], replace=FALSE),]
    offshore.sample<-offshore[sample(1:nrow(offshore), df12[i,2], replace=FALSE),]
    sample<-rbind(fringe.sample, offshore.sample)
    fish<-merge(sample, MNFishSeqRelSeqBPnc, by=c("SiteID","Depth","Habitat"))
    fish.2<-summaryBy(Presence.x~Species.Common.Name.x, data=fish, FUN=c(sum))
    fish.2$Presence_Absence.sum<-ifelse(fish.2$Presence.x.sum>=1,1,0)
    vec1[j]<-sum(fish.2$Presence_Absence.sum)
  }
  df12[i,4]<-mean(vec1)
}
#GRAPH OF RANDOMIZATION RESULTS
df12$prop.fringe<-df12$n.fringe/20
df12$prop.offshore<-df12$n.offshore/20
ggplot(df12, aes(prop.fringe, prop.offshore))+
  geom_tile(aes(fill = result), color='white')+
  scale_fill_gradient2(name="", low ="#0072B2", mid='white', high ="red4", space = 'rgb', guide = "colourbar", midpoint = (max(df12$result)+min(df12$result))/2, breaks=c(min(df12$result),(max(df12$result)+min(df12$result))/2, max(df12$result)), labels=c(round(min(df12$result)), round((max(df12$result)+min(df12$result))/2), round(max(df12$result))))+ 
  annotate("text", label="Total", x=ifelse(round(max(df12$result))>=10,0.93,0.95), y=1.0, hjust="right", colour="black")+ 
  annotate("rect", xmin=0.3, xmax=0.4, ymin=0.3, ymax=0.4, alpha=0.5)+
  theme_bw()+
  theme(legend.justification=c(1,1), legend.position=c(1,1), legend.title=element_blank())+
  scale_x_continuous(breaks=c(0,0.2,0.4, 0.6, 0.8, 1))+
  scale_y_continuous(breaks=c(0,0.2,0.4, 0.6, 0.8, 1))










#Move to DNR Dataset
#Limit DNR dataset to most recent population assessment that used all gear types
DNRFishRecentPA<-subset(DNRFish,Surveys=="Loon  06/05/2017 GN SD" |
                        Surveys=="Loon  06/05/2017 TN SD" |
                          Surveys=="Little Pine  06/28/2010 EFB PA" |
                          Surveys=="Little Pine  06/28/2010 GN PA" |
                          Surveys=="Little Pine  06/28/2010 S18 PA" |
                          Surveys=="Little Pine  06/28/2010 TN PA" |
                          Surveys=="Ely  08/07/2017 GN SD" |
                          Surveys=="Ely  08/07/2017 TN SD" |
                          Surveys=="East Chub  07/05/2022 TN SD" |
                          Surveys=="East Chub  07/05/2022 GN SD" |
                          Surveys=="Dunnigan  08/02/2021 GN SD" |
                          Surveys=="Dunnigan  08/02/2021 TN SD" |
                          Surveys=="Clark  06/01/2021 TN SD" |
                          Surveys=="Clark  06/01/2021 GN SD" |
                          Surveys=="Big Pine  08/30/2021 TN SD" |
                          Surveys=="Big Pine  08/30/2021 GN SD" |
                          Surveys=="Big Pine  07/06/2021 S58 TS" |
                          Surveys=="Big Pine  07/06/2021 EFB TS" |
                          Surveys=="Big Pine  08/30/2021 GN SD")

#Create species by site matrix
DNRPresMatrix<- cast(DNRFishRecentPA, Survey.ID.Date + Lake.Name + Station.Abbreviation...Number ~ Species.Common.Name, value='Presence')
#replace NA's with 0's
DNRPresMatrix[is.na(DNRPresMatrix)] <- 0
#melt back so there are 0's
DNRPres0s<-melt(DNRPresMatrix, na.rm = FALSE, value.name ="Presence", id =c("Survey.ID.Date","Lake.Name","Station.Abbreviation...Number"))

#aggregate to average number of samples required to get signal per lake
aggmeanDNRPres0s<- aggregate(DNRPres0s$value, by=list(DNRPres0s$Lake.Name,DNRPres0s$Species.Common.Name), FUN=mean)
#replace 0 with NA for visualization
aggmeanDNRPres<-subset(aggmeanDNRPres0s,x>0)
#visualize with heatmap
ggplot(aggmeanDNRPres, aes(Group.1, Group.2))+
  geom_tile(aes(fill = x))+
  scale_fill_gradient(low = "#ffffd9", high = "#081d58", name="Mean\n# Positive\nSamples")+
  labs(x="Lake", y="Common Name")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))

#pick species for AFS talk
#Walleye = game fish low abundance, DNR better
#Bluegill = game fish high abundance consistent
#pumpkinseed = similar but detected low abundance
#northern pike = game fish varied abundance, similar but detected low abundance
#largemouth bass = game fish same
#Round goby = invasive 
#Ruffe = invasive
#tubenose goby = invasive
#banded killifish = baitfish
#blacknose shiner = baitfish
#bluntnose minnow = baitfish
#brook stickleback = baitfish, very different detection

#Dig deeper into walleye
#Create DNR matrix with abundance rather than presence
#Create species by site matrix
DNRAbundMatrix<- cast(DNRFishRecentPA, Survey.ID.Date + Lake.Name + Station.Abbreviation...Number + Station.Type.Abbreviation ~ Species.Common.Name, value='Station.Total.Fish.Count')
#replace NA's with 0's
DNRAbundMatrix[is.na(DNRAbundMatrix)] <- 0
#melt back so there are 0's
DNRAbund0s<-melt(DNRAbundMatrix, na.rm = FALSE, value.name ="Abundance", id =c("Survey.ID.Date","Lake.Name","Station.Abbreviation...Number"))

#Create mean walleye abundance DNR
DNRAbund0sWalleye<-subset(DNRAbund0s, Species.Common.Name=="walleye")
#create mean lake numbers for walleye DNR for different gear methods
DNRAbund0sWalleyesummary<-summarySE(DNRAbund0sWalleye, measurevar="value",  groupvars=c("Lake.Name","Station.Type.Abbreviation","Survey.ID.Date"))
DNRAbund0sWalleyesummaryNets<-subset(DNRAbund0sWalleyesummary,Station.Type.Abbreviation=="GN"| 
                                       Station.Type.Abbreviation=="TN")
ggplot(DNRAbund0sWalleyesummaryNets, aes(x=Lake.Name, y=value,
                                         colour=Station.Type.Abbreviation)) +
  geom_point(size=2,position=position_dodge(width=0.8)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.4,position=position_dodge(width=0.8)) +  
  xlab("Lake")+
  ylab("CPUE")+
  scale_color_manual(values=c("black","orange"),name="Gear",
                    labels=c("Gill Net","Trap Net"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),legend.title=element_text(size=14))

#Create mean relative abundance DNA
MNFishSeqrelMatrix<- cast(MNFishSeqRelSeq, Lake.Name + SiteID + Depth + Habitat  ~ Species.Common.Name, value='relAbund')
#replace NA's with 0's
MNFishSeqrelMatrix[is.na(MNFishSeqrelMatrix)] <- 0
#melt back so there are 0's
MNFishSeqrel0s<-melt(MNFishSeqrelMatrix, na.rm = FALSE, value.name ="relAbund", id =c("Lake.Name","SiteID","Depth","Habitat"))

#Create mean walleye relative abundance DNA
MNFishSeqrel0sWalleye<-subset(MNFishSeqrel0s, Species.Common.Name=="Walleye" & Habitat!="choice")
#create mean lake numbers for walleye DNR for different gear methods
MNFishSeqrel0sWalleyesummary<-summarySE(MNFishSeqrel0sWalleye, measurevar="value",  groupvars=c("Lake.Name","Habitat","Depth"))
ggplot(MNFishSeqrel0sWalleyesummary, aes(x=Lake.Name, y=value,colour=Habitat)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.4) +  
  xlab("Lake")+
  ylab("Walleye Relative Sequence Abundance")+
  scale_colour_manual(name = "Sample Location", values = eDNASampleTypesNoChoiceColVec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.text=element_text(size=12),legend.title=element_blank())

#Create mean abundance DNA
MNFishSeqMatrix<- cast(MNFishSeqRelSeq, Lake.Name + SiteID + Depth + Habitat  ~ Species.Common.Name, value='SeqAbundance')
#replace NA's with 0's
MNFishSeqMatrix[is.na(MNFishSeqMatrix)] <- 0
#melt back so there are 0's
MNFishSeq0s<-melt(MNFishSeqMatrix, na.rm = FALSE, value.name ="SeqAbundance", id =c("Lake.Name","SiteID","Depth","Habitat"))

#Create mean walleye relative abundance DNA
MNFishSeq0sWalleye<-subset(MNFishSeq0s, Species.Common.Name=="Walleye" & Habitat!="choice")
#create mean lake numbers for walleye DNR for different gear methods
MNFishSeq0sWalleyesummary<-summarySE(MNFishSeq0sWalleye, measurevar="value",  groupvars=c("Lake.Name","Habitat","Depth"))
ggplot(MNFishSeq0sWalleyesummary, aes(x=Lake.Name, y=value,colour=Habitat)) +
  geom_point(size=2,position=position_dodge(width=0.8)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.5,position=position_dodge(width=0.8)) +  
  xlab("Lake")+
  ylab("Walleye Sequence Abundance")+
  scale_colour_manual(name = "Habitat", values =eDNASampleTypesNoChoiceColVec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.text=element_text(size=12),legend.title=element_text(size=14))

#correlation between sum of DNR and DNA
DNRAbu0sWalSumTotGN<-summarySE(subset(DNRAbund0sWalleye,Station.Type.Abbreviation=="GN"), measurevar="value",  groupvars=c("Lake.Name"))
DNRAbu0sWalSumTotGN<-DNRAbu0sWalSumTotGN%>% 
  rename(DNRN=N, DNRMean=value,DNRsd=sd,DNRse=se,DNRci=ci)
MNFishSeq0sWalSumTot<-summarySE(MNFishSeq0sWalleye, measurevar="value",  groupvars="Lake.Name")
MNFishSeq0sWalSumTot<-MNFishSeq0sWalSumTot%>% 
  rename(DNAN=N, DNAMean=value,DNAsd=sd,DNAse=se,DNAci=ci)
WalleyeCPUE<-merge(DNRAbu0sWalSumTotGN,MNFishSeq0sWalSumTot)
WalleyeCPUE$LakeClass<-WalleyeCPUE$Lake.Name
levels(WalleyeCPUE$LakeClass) <- list("Small"="East Chub", 
                                      "Small"="Dunnigan",
                                      "Medium"="Little Pine",
                                      "Medium"="Loon",
                                      "Medium"="Clark",
                                      "Large"="Big Pine",
                                      "Large"="Ripple",
                                      "Large"="Ely")
ggplot(WalleyeCPUE, aes(x=DNAMean, y=DNRMean,colour=LakeClass)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=DNRMean-DNRse, ymax=DNRMean+DNRse)) +  
  geom_errorbarh(aes(xmin=DNAMean-DNAse, xmax=DNAMean+DNAse)) +
  labs(x="eDNA CPUE", y="Gill Net CPUE")+
  scale_color_manual(name="Lake\nClass",values=c(LakeClassColVec))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.text=element_text(size=12),legend.title=element_text(size=14))

#Northern Pike
#Create mean walleye abundance DNR
DNRAbund0sNP<-subset(DNRAbund0s, Species.Common.Name=="northern pike")
#create mean lake numbers for walleye DNR for different gear methods
DNRAbund0sNPsummary<-summarySE(DNRAbund0sNP, measurevar="value",  groupvars=c("Lake.Name","Station.Type.Abbreviation","Survey.ID.Date"))
DNRAbund0sNPsummaryNets<-subset(DNRAbund0sNPsummary,Station.Type.Abbreviation=="GN"| 
                                       Station.Type.Abbreviation=="TN")
ggplot(DNRAbund0sNPsummaryNets, aes(x=Lake.Name, y=value,
                                         colour=Station.Type.Abbreviation)) +
  geom_point(size=2,position=position_dodge(width=0.8)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.4,position=position_dodge(width=0.8)) +  
  xlab("Lake")+
  ylab("CPUE")+
  scale_color_manual(values=c("black","orange"),name="Gear",
                     labels=c("Gill Net","Trap Net"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=12),
        legend.text=element_text(size=12),legend.title=element_text(size=14))

#Create mean walleye relative abundance DNA
MNFishSeqrel0sNP<-subset(MNFishSeqrel0s, Species.Common.Name=="Northern Pike" & Habitat!="choice")
#create mean lake numbers for walleye DNR for different gear methods
MNFishSeqrel0sNPsummary<-summarySE(MNFishSeqrel0sNP, measurevar="value",  groupvars=c("Lake.Name","Habitat","Depth"))
ggplot(MNFishSeqrel0sNPsummary, aes(x=Lake.Name, y=value,colour=Habitat)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.4) +  
  xlab("Lake")+
  ylab("Northern Pike Relative Sequence Abundance")+
  scale_colour_manual(name = "Sample Location", values = eDNASampleTypesNoChoiceColVec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.text=element_text(size=12),legend.title=element_blank())+
  facet_grid(Depth~.)

#Create mean northern pike relative abundance DNA
MNFishSeq0sNP<-subset(MNFishSeq0s, Species.Common.Name=="Northern Pike" & Habitat!="choice")
#create mean lake numbers for walleye DNR for different gear methods
MNFishSeq0sNPsummary<-summarySE(MNFishSeq0sNP, measurevar="value",  groupvars=c("Lake.Name","Habitat","Depth"))
ggplot(MNFishSeq0sNPsummary, aes(x=Lake.Name, y=value,colour=Habitat)) +
  geom_point(size=2,position=position_dodge(width=0.8)) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.5,position=position_dodge(width=0.8)) +  
  xlab("Lake")+
  ylab("Walleye Sequence Abundance")+
  scale_colour_manual(name = "Habitat", values =eDNASampleTypesNoChoiceColVec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.text=element_text(size=12),legend.title=element_text(size=14))+
  facet_grid(Depth~.)

#correlation between sum of DNR and DNA
DNRAbu0sNPSumTotGN<-summarySE(subset(DNRAbund0sNP,Station.Type.Abbreviation=="GN"), measurevar="value",  groupvars=c("Lake.Name"))
DNRAbu0sNPSumTotGN<-DNRAbu0sNPSumTotGN%>% 
  rename(DNRN=N, DNRMean=value,DNRsd=sd,DNRse=se,DNRci=ci)
MNFishSeq0sNPSumTot<-summarySE(MNFishSeq0sNP, measurevar="value",  groupvars="Lake.Name")
MNFishSeq0sNPSumTot<-MNFishSeq0sNPSumTot%>% 
  rename(DNAN=N, DNAMean=value,DNAsd=sd,DNAse=se,DNAci=ci)
NPCPUE<-merge(DNRAbu0sNPSumTotGN,MNFishSeq0sNPSumTot)
NPCPUE$LakeClass<-NPCPUE$Lake.Name
levels(NPCPUE$LakeClass) <- list("Small"="East Chub", 
                                      "Small"="Dunnigan",
                                      "Medium"="Little Pine",
                                      "Medium"="Loon",
                                      "Medium"="Clark",
                                      "Large"="Big Pine",
                                      "Large"="Ely")
ggplot(NPCPUE, aes(x=DNAMean, y=DNRMean,colour=LakeClass)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=DNRMean-DNRse, ymax=DNRMean+DNRse)) +  
  geom_errorbarh(aes(xmin=DNAMean-DNAse, xmax=DNAMean+DNAse)) +
  labs(x="eDNA CPUE", y="Gill Net CPUE")+
  scale_color_manual(name="Lake\nClass",values=c(LakeClassColVec))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.text=element_text(size=12),legend.title=element_text(size=14))


