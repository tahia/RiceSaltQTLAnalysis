library(qtl)
library(qtl2)
library(qtl2convert)
library(devtools)
library(devtools)
setwd("/home/taslima/data/JuengerLab/github/My_Gits/RiceSaltQTLAnalysis") #Choose directory 

rm(list=ls())
load("Reproductive_Reload_V5_RQTL2.RData")

set.seed(1234)

cross$pheno$SEL<-ifelse(cross$pheno$SES <=5, 1,ifelse(cross$pheno$SES >= 7, 3,2))
cross$pheno$Grp<-paste(cross$pheno$Cyto_num,cross$pheno$SEL,sep = "")

cross2<-convert2cross2(cross)
map <- insert_pseudomarkers(cross2, step=1)
pr <- calc_genoprob(cross2, map, error_prob=0.002, cores=16)


############## TEST FOR SELECTION ONLY FOR REPRODUCTIVE STAGE
covar <- as.data.frame(cross2$pheno[,c(2,13)])
strata<-(covar$Cyto_num>1)
#strata<-cross$pheno$Grp
names(strata)<-cross$pheno$id

kinship_loco <- calc_kinship(pr, "loco")

#Interactive
sc1_rep_loco_tail_int <- scan1(pr, cross2$pheno, kinship_loco, cores=16,addcovar = covar,intcovar = covar)
persc1_rep_loco_tail_int <- scan1perm(pr, cross2$pheno, kinship_loco ,n_perm=1000, cores=16,addcovar = covar,perm_strata = strata,intcovar = covar)

#Additive
sc1_rep_loco_tail_add <- scan1(pr, cross2$pheno, kinship_loco, cores=16,addcovar = covar)
persc1_rep_loco_tail_add <- scan1perm(pr, cross2$pheno, kinship_loco ,n_perm=1000, cores=16,addcovar = covar,perm_strata = strata)

#Null
covar <- as.data.frame(cross2$pheno[,c(2)])
names(covar)<-"Cyto_num"
strata<-(covar$Cyto_num>1)
names(strata)<-cross2$pheno[,1]

sc1_rep_loco_tail_null <- scan1(pr, cross2$pheno, kinship_loco, cores=16,addcovar = covar)
persc1_rep_loco_tail_null <- scan1perm(pr, cross2$pheno, kinship_loco ,n_perm=1000, cores=16,addcovar = covar,perm_strata = strata)

##### get threshold for selection
sel_threshold<-c()

for(i in 3:11)
{
  sel_threshold<-append(sel_threshold, as.numeric(quantile((persc1_rep_loco_tail_int[,i]- persc1_rep_loco_tail_null[,i]),probs = 0.95)))
}

###### find any QTL peak that pass this threshold
sel_sc1_rep_loco<-sc1_rep_loco_tail_int- sc1_rep_loco_tail_null
find_peaks(sel_sc1_rep_loco[,c(3:11)],map,sel_threshold)
find_peaks(sc1_rep_loco_tail_null,map,threshold = as.vector(summary(persc1_rep_loco_tail_null,alpha=0.05)))

###### Only Plant height 5@154 shows the SEL effect
############################### END OF TEST FOR THE EFFECT OF SELECTION ON REPRODUCTIVE STAGE

###### Now Start testing for additive and interactive effect of Cytoplasm as covariate

covar <- as.data.frame(cross2$pheno[,c(2)])
names(covar)<-"Cyto_num"
strata<-(covar$Cyto_num>1)
names(strata)<-cross2$pheno[,1]

#Interactive
sc1_rep_loco_cyt_int <- scan1(pr, cross2$pheno, kinship_loco, cores=16,addcovar = covar,intcovar = covar)
persc1_rep_loco_cyt_int <- scan1perm(pr, cross2$pheno, kinship_loco ,n_perm=1000, cores=16,addcovar = covar,perm_strata = strata,intcovar = covar)

#Additive
sc1_rep_loco_cyt_add <- scan1(pr, cross2$pheno, kinship_loco, cores=16,addcovar = covar)
persc1_rep_loco_cyt_add <- scan1perm(pr, cross2$pheno, kinship_loco ,n_perm=1000, cores=16,addcovar = covar,perm_strata = strata)

#Null
sc1_rep_loco_cyt_null <- scan1(pr, cross2$pheno, kinship_loco, cores=16)
persc1_rep_loco_cyt_null <- scan1perm(pr, cross2$pheno, kinship_loco ,n_perm=1000, cores=16)


##### get threshold for Interaction model where Cytoplasm is covar
int_cyt_threshold<-c()

for(i in 3:11)
{
  int_cyt_threshold<-append(int_cyt_threshold, as.numeric(quantile((persc1_rep_loco_cyt_int[,i]- persc1_rep_loco_cyt_null[,i]),probs = 0.95)))
}

###### find any QTL peak that pass this threshold for Cyt int
cyt_int_sc1_rep_loco<-sc1_rep_loco_cyt_int- sc1_rep_loco_cyt_null
find_peaks(cyt_int_sc1_rep_loco[,c(3:11)],map,int_cyt_threshold)


########### PH QTL at chr5@172, FGN, FGW at chr10@58.48, SF at chr10@60, HI at chr 10@59 shows interaction that cross the threshold
#### These positions are relative to models there will move a bit up and down


##### get threshold for Additive model where Cytoplasm is covar
add_cyt_threshold<-c()

for(i in 3:11)
{
  add_cyt_threshold<-append(add_cyt_threshold, as.numeric(quantile((persc1_rep_loco_cyt_add[,i]- persc1_rep_loco_cyt_null[,i]),probs = 0.95)))
}

###### find any QTL peak that pass this threshold for Cyt int
cyt_add_sc1_rep_loco<-sc1_rep_loco_cyt_add- sc1_rep_loco_cyt_null
find_peaks(cyt_add_sc1_rep_loco[,c(3:11)],map,add_cyt_threshold)

####### HI QTL at chr10@105 shows additive effect that cross the threshold

##So, overall for all all the traits for reproductive phenotype Cytoplasm will be 
##considered as both additive and interactive cofactor to test for QTl
## while asking for effect of a single QTL for a given position linear model will
##be fitted and Cytoplasm as covariates will be included if and only if it pass the 
##threshold of the difference between full and reduced model

###########################

RepQ<-as.data.frame(find_peaks(sc1_rep_loco_cyt_int,map,as.vector(summary(persc1_rep_loco_cyt_int,alpha=0.05))))

color <- c("slateblue", "violetred", "green3")
for(i in 3:11) {
  ymx <- max(c(as.data.frame(sc1_rep_loco_cyt_int)[,i],as.data.frame(sc1_rep_loco_cyt_add)[,i],as.data.frame(sc1_rep_loco_cyt_null)[,i]),na.rm =T) # overall maximum LOD score
  plot(sc1_rep_loco_cyt_null, map, lodcolumn=i, col=color[1], main=colnames(cross2$pheno)[i],
       ylim=c(0, ymx*1.02))
  plot(sc1_rep_loco_cyt_add, map, lodcolumn=i, col=color[2], add=TRUE)
  plot(sc1_rep_loco_cyt_int, map, lodcolumn=i, col=color[3], add=TRUE, lty=2)
  abline(h = as.vector(summary(persc1_rep_loco_cyt_null, alpha= 0.05))[i],col=color[1])
  abline(h = as.vector(summary(persc1_rep_loco_cyt_add, alpha= 0.05))[i],col=color[2])
  abline(h = as.vector(summary(persc1_rep_loco_cyt_int, alpha= 0.05))[i],col=color[3])
  legend("topleft", lwd=2, col=color, c("Null", "Cyt-Add", "Cyt-Int"), bg="gray90", lty=c(1,1,2))
}

save.image("Rep_QTL_Fin_v1.RData")

########## Now fit single QTL model
###### PH
covar_factor <- factor(covar$Cyto_num)
covar_matrix <- as.matrix(model.matrix( ~ covar_factor)[ , -1])

rep_PH_Q1<-fit1(pr[[1]][,,215],cross2$pheno[,'PH'],addcovar = covar,intcovar = covar)

test<-as.data.frame(attr(scan1coef(pr[,1],cross2$pheno[,'PH'],addcovar = covar,intcovar = covar,se = T),"SE"))
test<-scan1coef(pr[,1],cross2$pheno[,'PH'],addcovar = covar,intcovar = covar,se = T)

color <- c("slateblue", "violetred", "green3")
for(i in 3:11) {
  plot(test, map[1], columns =c(1:3), col=color, main=colnames(cross2$pheno)[i])
  plot(test$AB, map, lodcolumn=i, col=color[2], add=TRUE)
  plot(test$BB, map, lodcolumn=i, col=color[3], add=TRUE, lty=2)
  legend("topleft", lwd=2, col=color, c("AA", "AB", "BB"), bg="gray90", lty=c(1,1,2))
}

test2<-attr(scan1coef(pr[,1],cross2$pheno[,'PH'],se = T),"SE")


############ Seedling Stage
rm(list=ls())
set.seed(1234)
load("Seedling_Reload_V5_RQTL2.RData")

cross2<-convert2cross2(cross)
map <- insert_pseudomarkers(cross2, step=1)
pr <- calc_genoprob(cross2, map, error_prob=0.002, cores=16)
covar <- as.data.frame(cross2$pheno[,c(2)])
names(covar)<-"Cyto_num"
strata<-(covar$Cyto_num>1)
names(strata)<-cross2$pheno[,1]
kinship_loco <- calc_kinship(pr, "loco")

###### Now Start testing for additive and interactive effect of Cytoplasm as covariate

#Interactive
sc1_seed_loco_cyt_int <- scan1(pr, cross2$pheno, kinship_loco, cores=16,addcovar = covar,intcovar = covar)
persc1_seed_loco_cyt_int <- scan1perm(pr, cross2$pheno, kinship_loco ,n_perm=1000, cores=16,addcovar = covar,perm_strata = strata,intcovar = covar)

#Additive
sc1_seed_loco_cyt_add <- scan1(pr, cross2$pheno, kinship_loco, cores=16,addcovar = covar)
persc1_seed_loco_cyt_add <- scan1perm(pr, cross2$pheno, kinship_loco ,n_perm=1000, cores=16,addcovar = covar,perm_strata = strata)

#Null
sc1_seed_loco_cyt_null <- scan1(pr, cross2$pheno, kinship_loco, cores=16)
persc1_seed_loco_cyt_null <- scan1perm(pr, cross2$pheno, kinship_loco ,n_perm=1000, cores=16)


##### get threshold for Interaction model where Cytoplasm is covar
int_cyt_seed_threshold<-c()

for(i in 3:10)
{
  int_cyt_seed_threshold<-append(int_cyt_seed_threshold, as.numeric(quantile((persc1_seed_loco_cyt_int[,i]- persc1_seed_loco_cyt_null[,i]),probs = 0.95)))
}

###### find any QTL peak that pass this threshold for Cyt int
cyt_int_sc1_seed_loco<-sc1_seed_loco_cyt_int- sc1_seed_loco_cyt_null
find_peaks(cyt_int_sc1_seed_loco[,c(3:10)],map,int_cyt_seed_threshold)


########### 
#### These positions are relative to models there will move a bit up and down


##### get threshold for Additive model where Cytoplasm is covar
add_cyt_seed_threshold<-c()

for(i in 3:10)
{
  add_cyt_seed_threshold<-append(add_cyt_seed_threshold, as.numeric(quantile((persc1_seed_loco_cyt_add[,i]- persc1_seed_loco_cyt_null[,i]),probs = 0.95)))
}

###### find any QTL peak that pass this threshold for Cyt add
cyt_add_sc1_seed_loco<-sc1_seed_loco_cyt_add- sc1_seed_loco_cyt_null
find_peaks(cyt_add_sc1_seed_loco[,c(3:10)],map,add_cyt_seed_threshold)

#######  shows additive effect that cross the threshold

##So, overall for all all the traits for reproductive phenotype Cytoplasm will be 
##considered as both additive and interactive cofactor to test for QTl
## while asking for effect of a single QTL for a given position linear model will
##be fitted and Cytoplasm as covariates will be included if and only if it pass the 
##threshold of the difference between full and reduced model

###########################

SeedQ<-as.data.frame(find_peaks(sc1_seed_loco_cyt_int,map,as.vector(summary(persc1_seed_loco_cyt_int,alpha=0.05))))

color <- c("slateblue", "violetred", "green3")
for(i in 3:10) {
  ymx <- max(c(as.data.frame(sc1_seed_loco_cyt_int)[,i],as.data.frame(sc1_seed_loco_cyt_add)[,i],as.data.frame(sc1_seed_loco_cyt_null)[,i]),na.rm =T) # overall maximum LOD score
  plot(sc1_seed_loco_cyt_null, map, lodcolumn=i, col=color[1], main=colnames(cross2$pheno)[i],
       ylim=c(0, ymx*1.02))
  plot(sc1_seed_loco_cyt_add, map, lodcolumn=i, col=color[2], add=TRUE)
  plot(sc1_seed_loco_cyt_int, map, lodcolumn=i, col=color[3], add=TRUE, lty=2)
  abline(h = as.vector(summary(persc1_seed_loco_cyt_null, alpha= 0.05))[i],col=color[1])
  abline(h = as.vector(summary(persc1_seed_loco_cyt_add, alpha= 0.05))[i],col=color[2])
  abline(h = as.vector(summary(persc1_seed_loco_cyt_int, alpha= 0.05))[i],col=color[3])
  legend("topleft", lwd=2, col=color, c("Null", "Cyt-Add", "Cyt-Int"), bg="gray90", lty=c(1,1,2))
}

save.image("Seed_QTL_Fin_v1.RData")

###############################
#get flanking markers
dat<-read.csv("RiceSaltQTL_table2.csv",check.names = F)

colnames(dat)
#colnames(dat)[3]<-"Chromosome"
colnames(dat)[which(colnames(dat) == "Chr")]<-"Chromosome"
dat$'QTL Interval'<-paste(dat$`upstream position  of peak`,dat$`downstream position of peak`,sep = " - ")
dat<-dat[,c(1:5,9,8)]
colnames(dat)
colnames(dat)[which(colnames(dat) == "Position")]<-"QTL peak"


### Will use rqtl to get the real markers, not pseudo
# get_flanking_markers<-function(cross,dat){
#   up=down=c()
#   for (i in 1:nrow(dat)) {
#       p<-as.numeric(strsplit(dat[i,'QTL Interval'],split = " - ")[[1]]);
#       chr<-dat[i,'Chromosome']
#       library(qtl)
#       up<-append(up,strsplit(find.marker(cross,chr,p[1]),split = "_")[[1]][3])
#       down<-append(down,strsplit(find.marker(cross,chr,p[2]),split = "_")[[1]][3])
#       }
#       list=list(up=up,down=down)
#       return(list)
# }

############ Re-write the above function because some genetic vs physical order are not same


get_flanking_markers<-function(cross,dat){
  up=down=c()
  for (i in 1:nrow(dat)) {
    p<-as.numeric(strsplit(dat[i,'QTL Interval'],split = " - ")[[1]]);
    chr<-dat[i,'Chromosome']
    library(qtl)
    BP<-c()
    for (j in 1:length(unlist(pull.map(cross,chr)))) 
    { 
      BP<-append(BP,as.numeric(strsplit(names(unlist(pull.map(cross,chr))),split = "_")[[j]][3]) )
      }
    CM<-as.vector(unlist(pull.map(cross,chr)))
    genotab<-as.data.frame(cbind(BP=BP,CM=CM))
    up<-append(up,min(genotab$BP[which(genotab$CM >= round(p[1]) & genotab$CM < round(p[2]))]))
    down<-append(down,max(genotab$BP[which(genotab$CM >= round(p[1]) & genotab$CM < round(p[2]))]))
  }
  list=list(up=up,down=down)
  return(list)
}

fbps<-get_flanking_markers(cross,dat)

dat$upBP<-fbps$up
dat$downBP<-fbps$down

write.csv(dat,"RiceSaltQTL_table2_withflankingM.csv",row.names = F)

