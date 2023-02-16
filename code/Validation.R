library(tidyverse)
library(rstan)
library(loo)
library(rstanarm)
library(pROC)
library(survival)
library(survminer)
library(rms)

# METABRIC validation
#Read patient data downloaded from cbioportal 
patient_data <- read.csv("../data/brca_metabric/data_clinical_patient.txt")
colnames(patient_data) <- patient_data[1,]
patient_data <- patient_data[-1,]

HT_patients <- patient_data[which(patient_data$HORMONE_THERAPY=="YES"),]
HTERpos_patients <- HT_patients[which(HT_patients$ER_IHC=="Positve"),]
HT_CancerDeath_patients <- HTERpos_patients[-which(HTERpos_patients$VITAL_STATUS=="Died of Other Causes"),]

Neg_patients <- patient_data[which(patient_data$ER_IHC=="Negative"),]

# Read RNA-seq data downloaded from cbioportal
RNA_seq <- read.csv("../data/brca_metabric/data_mrna_agilent_microarray_zscores_ref_all_samples.txt", header = F, sep="\t")
colnames(RNA_seq) <- RNA_seq[1,]
RNA_seq <- data.frame(t(RNA_seq[-1,]))
colnames(RNA_seq) <- RNA_seq[1,]
RNA_seq <- RNA_seq[-c(1,2),]

#Adapt data to work with gene values
data_mrna_all <- RNA_seq
data_mrna_all$Patient <- rownames(data_mrna_all)
data_mrna_all[,-dim(data_mrna_all)[2]] <- data_mrna_all[-dim(data_mrna_all)[2]] %>% 
  apply( 2, as.numeric) %>% data.frame()
data_mrna_all <- data_mrna_all[,-duplicated(colnames(data_mrna_all))]


#Create KM dataframe for HT patients
km_data_all <- inner_join(data_mrna_all, HT_CancerDeath_patients, by=c("Patient"="PATIENT_ID"))

#Create KM dataframe for ER Negative
km_data_control <- inner_join(data_mrna_all, Neg_patients, by=c("Patient"="PATIENT_ID"))

km_data_all$RFS_MONTHS <- as.numeric(km_data_all$RFS_MONTHS)
km_data_all$RFS_MONTHS<- as.numeric(km_data_all$RFS_MONTHS)
km_data_all$RFS_STATUS <- gsub(x=km_data_all$RFS_STATUS, "0:Not Recurred", 0)
km_data_all$RFS_STATUS <- gsub(x=km_data_all$RFS_STATUS, "1:Recurred", 1)
km_data_all$RFS_STATUS <- as.numeric(km_data_all$RFS_STATUS)
km_data_all$RFS_STATUS <- ifelse(km_data_all$RFS_MONTHS > 120, 0, km_data_all$RFS_STATUS)
km_data_all$RFS_MONTHS <- ifelse(km_data_all$RFS_MONTHS > 120, 120, km_data_all$RFS_MONTHS)

km_data_control$RFS_MONTHS <- as.numeric(km_data_control$RFS_MONTHS)
km_data_control$RFS_MONTHS<- as.numeric(km_data_control$RFS_MONTHS)
km_data_control$RFS_STATUS <- gsub(x=km_data_control$RFS_STATUS, "0:Not Recurred", 0)
km_data_control$RFS_STATUS <- gsub(x=km_data_control$RFS_STATUS, "1:Recurred", 1)
km_data_control$RFS_STATUS <- as.numeric(km_data_control$RFS_STATUS)
km_data_control$RFS_STATUS <- ifelse(km_data_control$RFS_MONTHS > 120, 0, km_data_control$RFS_STATUS)
km_data_control$RFS_MONTHS <- ifelse(km_data_control$RFS_MONTHS > 120, 120, km_data_control$RFS_MONTHS)

#### Generation of Survival curves and Cox Models ####


##Our 6 gene signature

# ER positive, HT treated (non-cancer related death excluded) patients
OurSig <- c("VTN", "GRB14", "HERC1", "FRAS1", "INSIG2","TMC7")

km_data_all$Oursig <- apply(km_data_all[,OurSig],1,mean)

Upper_t_Ours <- as.numeric(quantile(km_data_all$Oursig, probs = seq(0, 1, 1/3))[3], type=1)

km_data_all$Risk_Ours <- ifelse(km_data_all$Oursig > Upper_t_Ours,
                                "High","Medium/Low")

km_data_all$Risk_Ours <- factor(km_data_all$Risk_Ours,
                                levels = c("Medium/Low", "High"))

cph_Ours <- coxph(Surv(km_data_all$RFS_MONTHS,
                       km_data_all$RFS_STATUS) ~ Risk_Ours, km_data_all)
summary(cph_Ours)

fit_Ours <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ Risk_Ours,
                    km_data_all)

sd_Ours <- survdiff(Surv(km_data_all$RFS_MONTHS,
                         km_data_all$RFS_STATUS ) ~ Risk_Ours, km_data_all)
p_val_Ours =1 - pchisq(sd_Ours$chisq, length(sd_Ours$n) - 1)

p_val_Ours
summary(cph_Ours)


a <- survminer::ggsurvplot(
  fit = fit_Ours, 
  palette = c("black", "red"),
  xlab = "Time (months)", ylab = "RFS Probability",
  #title = "Metabric ER+ HT treated patients",
  risk.table=T,  risk.table.col = "strata",
  #pval,
  break.x.by=20, break.y.by = 0.2, 
  legend=c(0.13, 0.2),
  legend.title = "Expression",
  legend.labs = c("low", "high"),
  ggtheme=theme_linedraw())

a$table <- a$table + theme(axis.line = element_blank(), panel.grid = element_blank(), 
                           text = element_text(size = 15),axis.title = element_blank(),
                           axis.text = element_blank(), panel.border = element_blank()) 
a$plot <-  a$plot + theme(text = element_text(size = 20), panel.grid = element_blank()) +
  annotate("text", x=90, y=0.15, size=7,
           label= " HR = 1.57 (1.25 - 1.96)\n      logrank P = 6.7e-05") 


print(a)


## Negative control


# ER negative patients
OurSig <- c("VTN", "GRB14", "HERC1", "FRAS1", "INSIG2","TMC7")

km_data_control$Oursig <- apply(km_data_control[,OurSig],1,mean)

Upper_t_Ours <- as.numeric(quantile(km_data_control$Oursig, probs = seq(0, 1, 1/3))[3], type=1)

km_data_control$Risk_Ours <- ifelse(km_data_control$Oursig > Upper_t_Ours,
                                "High","Medium/Low")

km_data_control$Risk_Ours <- factor(km_data_control$Risk_Ours,
                                levels = c("Medium/Low", "High"))

cph_Ours <- coxph(Surv(km_data_control$RFS_MONTHS,
                       km_data_control$RFS_STATUS) ~ Risk_Ours, km_data_control)
summary(cph_Ours)

fit_Ours <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ Risk_Ours,
                    km_data_control)

sd_Ours <- survdiff(Surv(km_data_control$RFS_MONTHS,
                         km_data_control$RFS_STATUS ) ~ Risk_Ours, km_data_control)
p_val_Ours =1 - pchisq(sd_Ours$chisq, length(sd_Ours$n) - 1)

p_val_Ours
summary(cph_Ours)


a <- survminer::ggsurvplot(
  fit = fit_Ours, 
  palette = c("black", "red"),
  xlab = "Time (months)", ylab = "RFS Probability",
  #title = "Metabric ER+ HT treated patients",
  risk.table=T,  risk.table.col = "strata",
  #pval,
  break.x.by=20, break.y.by = 0.2, 
  legend=c(0.13, 0.2),
  legend.title = "Expression",
  legend.labs = c("low", "high"),
  ggtheme=theme_linedraw())

a$table <- a$table + theme(axis.line = element_blank(), panel.grid = element_blank(), 
                           text = element_text(size = 15),axis.title = element_blank(),
                           axis.text = element_blank(), panel.border = element_blank()) 
a$plot <-  a$plot + theme(text = element_text(size = 20), panel.grid = element_blank())


print(a)




## 5 candidate pathways
# For Rahem2020, get ROC-UAC values as they are used as weights in the signature
data_Raheem <- km_data_all

Get_ROC <- function(Gene, data_Raheem){
  data <- data.frame(data_Raheem[,Gene])
  data$Response <- data_Raheem$RFS_STATUS
  
  fit <- stan_glm(formula = Response ~ .,
                  data = data, family = binomial(link = "logit"),
                  refresh = 0)
  
  loo <- loo(fit, save_psis = T)
  preds <- posterior_epred(fit)
  pred <- colMeans(preds)
  ploo = E_loo(preds, loo$psis_object, type="mean",
               log_ratios = -log_lik(fit))$value
  
  ROC <- roc(data$Response,ploo, plot = F) #, legacy.axes = TRUE,ci=TRUE, #smooth=TRUE,
  #percent = F,  col = "#477eb9", lwd = 2, print.auc = TRUE)
  return(ROC)
}

LARS_weight <- Get_ROC("LARS", data_Raheem)
AP2S1_weight <-Get_ROC("AP2S1", data_Raheem)
GTF3C3_weight <- Get_ROC("GTF3C3", data_Raheem)
EIF2AK3_weight <- Get_ROC("EIF2AK3", data_Raheem)
CDK1_weight <- Get_ROC("CDK1", data_Raheem)


data_Raheem$Risk_Score <- (LARS_weight$auc*data_Raheem$LARS +
                             AP2S1_weight$auc*data_Raheem$AP2S1 +
                             GTF3C3_weight$auc*data_Raheem$GTF3C3 +
                             EIF2AK3_weight$auc*data_Raheem$EIF2AK3 +
                             CDK1_weight$auc*data_Raheem$CDK1)

#Create risk groups and make CPH calculation
mean_risk <- mean(data_Raheem$Risk_Score)
sd_risk <- sd(data_Raheem$Risk_Score)
data_Raheem$Group_risk <- ifelse(data_Raheem$Risk_Score > (mean_risk + sd_risk),
                                 "High Risk", "Medium/Low Risk")
data_Raheem$Risk_Raheem <- data_Raheem$Group_risk

cph_Raheem <- coxph(Surv(RFS_MONTHS, RFS_STATUS) ~ Risk_Raheem,
                    data=data_Raheem)
summary(cph_Raheem)

fit_Raheem <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ Risk_Raheem, 
                      data=data_Raheem)

sd_Raheem <- survdiff(Surv(RFS_MONTHS, RFS_STATUS) ~ Risk_Raheem,
                      data_Raheem)
p_val_Raheem =1 - pchisq(sd_Raheem$chisq, length(sd_Raheem$n) - 1)

ggsurvplot(
  fit = fit_Raheem, 
  palette = c("red", "black"),
  xlab = "Months", ylab = "Relapse Free Survival",
  title = "Metabric ER+", subtitle="Raheem et al. Gene Signature",  
  risk.table=T, tables.y.text=F, 
  pval=p_val_Raheem,break.time.by=24,
  ggtheme=theme_survminer(font.main = 20, font.x = 20, font.y = 20, 
                          font.submain = 15, font.tickslab = 20, 
                          font.legend = 20))

DXY_Raheem <- cph(Surv(RFS_MONTHS, RFS_STATUS) ~ Risk_Raheem,
                  data_Raheem, x=T, y=T)
validate(DXY_Raheem, method='boot', B=150)

############ SET ER PR 
T_genes <- c("SLC39A6", "STC2", "CA12", "ESR1", "PDZK1", "NPY1R", "CD2", "MAPT", "QDPR",
             "AZGP1", "ABAT", "ADCY1", "CD3D", "NAT1", "MRPS30", "DNAJC12", "SCUBE2", "KCNE4")

R_genes <- c("LDHA", "ATP5J2", "VDAC2", "DARS", "UGP2", "UBE2Z", "AK2", "WIPF2", "APPBP2", "TRIM2")

Mean_T <- apply(km_data_all[,T_genes],1,mean)
Mean_R <- apply(km_data_all[,R_genes],1,mean)

km_data_all$SET <- Mean_T - Mean_R + 2

km_data_all$Risk_SET <- ifelse(km_data_all$SET <  median(km_data_all$SET),
                               "High","Medium/Low")

cph_SET<- coxph(Surv(km_data_all$RFS_MONTHS,
                     km_data_all$RFS_STATUS) ~ Risk_SET, km_data_all)
summary(cph_SET)

fit_SET <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ Risk_SET,
                   km_data_all)

sd_SET <- survdiff(Surv(km_data_all$RFS_MONTHS,
                        km_data_all$RFS_STATUS ) ~ Risk_SET, km_data_all)
p_val_SET =1 - pchisq(sd_SET$chisq, length(sd_SET$n) - 1)

ggsurvplot(
  fit = fit_SET, 
  palette = c("red", "black"),
  xlab = "Months", ylab = "Relapse Free Survival",
  title = "Metabric ER+ HT treated patients", subtitle="GSS - Upper tertile division",  
  risk.table=T, tables.y.text=F, 
  pval=p_val_SET, break.time.by=24,
  ggtheme=theme_survminer(font.main = 20, font.x = 20, font.y = 20, 
                          font.submain = 15, font.tickslab = 20, 
                          font.legend = 20))

DXY_SET <- cph(Surv(RFS_MONTHS, RFS_STATUS) ~ Risk_SET,
               km_data_all, x=T, y=T)
validate(DXY_SET, method='boot', B=150)

####### CRISPR Sig
CRISPRSig <- c("MTMR7", "SMTNL2", "GATA4", "IL20RA", "CA12", "EPB49")


km_data_all$CrisprSig <- apply(km_data_all[,CRISPRSig],1,mean)

# Calculate optimal cut value for Low/High risk patients
Crispr_Cut<- surv_cutpoint(data =km_data_all, time = "RFS_MONTHS",
                           event = "RFS_STATUS", 
                           variables="CrisprSig"
)

Crispr_Risk <- surv_categorize(Crispr_Cut)


cph_Crispr <- coxph(Surv(RFS_MONTHS, RFS_STATUS) ~ CrisprSig,
                    data=Crispr_Risk )
summary(cph_Crispr)

fit_Crispr <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ CrisprSig, 
                      data=Crispr_Risk )

sd_Crispr <- survdiff(Surv(RFS_MONTHS, RFS_STATUS) ~ CrisprSig,
                      Crispr_Risk )
p_val_Crispr =1 - pchisq(sd_Crispr$chisq, length(sd_Crispr$n) - 1)

ggsurvplot(
  fit = fit_Crispr, 
  palette = c("red", "black"),
  xlab = "Months", ylab = "Relapse Free Survival",
  title = "Metabric ER+ HT treated patients", subtitle="CRISPR Gene Signature",  
  risk.table=T, tables.y.text=F, 
  pval=p_val_Crispr,break.time.by=24,
  ggtheme=theme_survminer(font.main = 20, font.x = 20, font.y = 20, 
                          font.submain = 15, font.tickslab = 20, 
                          font.legend = 20))


####### ODX score
ODXSig <- c("ACTB", "GAPDH","RPLP0","GUSB", "TFRC",
            "ESR1", "PGR", "BCL2", "SCUBE2","GRB7",
            "ERBB2","MKI67", "AURKA","BIRC5", "CCNB1",
            "MYBL2","MMP11","CTSL2", "GSTM1", "CD68")

#Taking the step previous to the join and joining later due to an error returning a blank dataset for Model_ODX
ODX_patients <- HT_CancerDeath_patients[which(HT_CancerDeath_patients$HER2_SNP6!="GAIN"),]

ODX_RS <- function (mat) {
  # Define list of genes
  Ref_grp <- c("ACTB", "GAPDH", "RPLP0", "GUSB", "TFRC")
  ER_grp <- c("ESR1", "PGR", "BCL2", "SCUBE2")
  HER2_grp <- c("GRB7", "ERBB2")
  Prol_grp <- c("MKI67", "AURKA", "BIRC5", "CCNB1", "MYBL2")
  Inv_grp <- c("MMP11", "CTSL2")
  
  All_grp <- c(Ref_grp, ER_grp, HER2_grp, Prol_grp, Inv_grp, "GSTM1", "CD68") # Missing BAG1 in the dataset
  
  # Filter for 21 genes  
  mat <- data.frame(mat[, match(All_grp, colnames(mat))])
  
  # Rescale expression values to lie in the range of 1:20
  mat <- (((mat - min(mat)) / (max(mat) - min(mat))) * 30) + 20
  
  # Normalize by reference genes 
  mat <- mat - rowMeans(mat[, Ref_grp])
  
  # Calculate group scores 
  Her2_score <- 0.9*mat$"GRB7" + 0.1*mat$"ERBB2"
  Her2_score[which(Her2_score < 8)]  <- 8 
  
  ER_score <- (0.8*mat$"ESR1" + 1.2*mat$"PGR" + mat$"BCL2" + mat$"SCUBE2") / 4
  
  Prol_score <- rowMeans(mat[ ,Prol_grp])
  Prol_score[which(Prol_score < 6.5)] <- 6.5
  
  Inv_score <- rowMeans(mat[, Inv_grp])
  
  # Calculate unscaled recurrence scores
  RSu <- 0.47*Her2_score + -0.34*ER_score + 1.04*Prol_score + 0.1*Inv_score + 0.05*mat$"CD68" + -0.08*mat$"GSTM1" #+ -0.07*mat$"BAG1"
  
  # Calculate scaled recurrence scores
  # if (RSu > 100) RS <- 100
  #  if (RSu < 0) RS <- 0
  #if (RSu > 0 && RSu < 100) 
  RS <- 20*(RSu - 6.7)
  return(RS)
}
ODX_genes <- data_mrna_all[,which(colnames(data_mrna_all) %in% ODXSig)]

ODX <- ODX_RS(ODX_genes)
ODX <- ifelse(ODX>100, 100,ODX)

Model_ODX <- data.frame(ODX = ODX, Patient =rownames(ODX_genes))

Model_ODX <- inner_join(Model_ODX, HT_CancerDeath_patients, by=c("Patient"="PATIENT_ID"))

Model_ODX$RFS_MONTHS <- as.numeric(Model_ODX$RFS_MONTHS)
Model_ODX$RFS_MONTHS<- as.numeric(Model_ODX$RFS_MONTHS)
Model_ODX$RFS_STATUS <- gsub(x=Model_ODX$RFS_STATUS, "0:Not Recurred", 0)
Model_ODX$RFS_STATUS <- gsub(x=Model_ODX$RFS_STATUS, "1:Recurred", 1)
Model_ODX$RFS_STATUS <- as.numeric(Model_ODX$RFS_STATUS)
Model_ODX$RFS_STATUS <- ifelse(Model_ODX$RFS_MONTHS > 120, 0, Model_ODX$RFS_STATUS)
Model_ODX$RFS_MONTHS <- ifelse(Model_ODX$RFS_MONTHS > 120, 120, Model_ODX$RFS_MONTHS)

# The 85% percentile for the division comes from the distribution in:
# Cheng, R., Kong, X., Wang, X., Fang, Y., & Wang, J. (2020). Oncotype DX breast recurrence score distribution 
# and chemotherapy benefit among women of different age groups with HR-positive,
# HER2-negative, node-negative breast cancer in the SEER database. 
# Frontiers in Oncology, 10, 1583.

#### 1389/9719 = 0.1429159% from:


Upper_t_ODX <- as.numeric(quantile(Model_ODX$ODX, probs = seq(0, 1, 1/5))[5], type=1)

Model_ODX$Group_risk <- ifelse(Model_ODX$ODX > Upper_t_ODX,
                               "High Risk", "Medium/Low Risk")
Model_ODX$Risk_ODX <- Model_ODX$Group_risk

cph_ODX <- coxph(Surv(RFS_MONTHS, RFS_STATUS) ~ ODX,
                 data=Model_ODX)
summary(cph_ODX)

fit_ODX <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ Risk_ODX, 
                   data=Model_ODX)

sd_ODX <- survdiff(Surv(RFS_MONTHS, RFS_STATUS) ~ Risk_ODX,
                   Model_ODX)
p_val_ODX =1 - pchisq(sd_ODX$chisq, length(sd_ODX$n) - 1)

ggsurvplot(
  fit = fit_ODX, 
  palette = c("red", "black"),
  xlab = "Months", ylab = "Relapse Free Survival",
  title = "Metabric ER+ HT patients", subtitle="ODX Score",  
  risk.table=T, tables.y.text=F, 
  pval=p_val_ODX,break.time.by=24,
  ggtheme=theme_survminer(font.main = 20, font.x = 20, font.y = 20, 
                          font.submain = 15, font.tickslab = 20, 
                          font.legend = 20))

DXY_ODX<- cph(Surv(RFS_MONTHS, RFS_STATUS) ~ Risk_ODX,
              Model_ODX, x=T, y=T)
validate(DXY_ODX, method='boot', B=150)

######## Men et al
MenSig <- c("GRB14","NPY1R","SOX8","GABRP","GFRA3",
            "PTPRN2","ARNT2","ATP6V1C2","NELL2") #,"C2CD4D" Not in METABRIC

km_data_all$Mensig <- apply(km_data_all[,MenSig],1,mean)

#USing 1/3 threshold 
# Upper_t_Mens <- as.numeric(quantile(km_data_all$Mensig, probs = seq(0, 1, 1/3))[3], type=1)
# km_data_all$Risk_Men <- ifelse(km_data_all$Mensig >= Upper_t_Mens,
#                                 "High","Medium/Low")


# Calculate optimal cut value for Low/High risk patients
Men_Cut<- surv_cutpoint(data =km_data_all, time = "RFS_MONTHS",
                        event = "RFS_STATUS", 
                        variables="Mensig"
)

Men_Risk <- surv_categorize(Men_Cut)


cph_Men <- coxph(Surv(RFS_MONTHS, RFS_STATUS) ~ Mensig,
                 data=Men_Risk )

# cph_Men<- coxph(Surv(km_data_all$RFS_MONTHS,
#                      km_data_all$RFS_STATUS) ~ Risk_Men, km_data_all)
# summary(cph_Men)

fit_Men <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ Mensig,
                   Men_Risk)

sd_Men <- survdiff(Surv(km_data_all$RFS_MONTHS,
                        km_data_all$RFS_STATUS ) ~ Mensig, Men_Risk)
p_val_Men =1 - pchisq(sd_Men$chisq, length(sd_Men$n) - 1)

ggsurvplot(
  fit = fit_Men, 
  palette = c("red", "black"),
  xlab = "Months", ylab = "Relapse Free Survival",
  title = "Metabric ER+ HT treated patients", subtitle="Men et al - 10 gene signature",  
  risk.table=T, tables.y.text=F, 
  pval=p_val_Men, break.time.by=24,
  ggtheme=theme_survminer(font.main = 20, font.x = 20, font.y = 20, 
                          font.submain = 15, font.tickslab = 20, 
                          font.legend = 20))

DXY_Men<- cph(Surv(RFS_MONTHS, RFS_STATUS) ~ Risk_Men,
              km_data_all, x=T, y=T)
validate(DXY_Men, method='boot', B=150)

######## Ma et al
MaSig <- c("HOXB13","IL17RB")

km_data_all$Risk_Ma <- ifelse( km_data_all$HOXB13 > km_data_all$IL17RB,
                               "High","Medium/Low")
km_data_all$Risk_Ma <- as.factor(km_data_all$Risk_Ma)

cph_Ma<- coxph(Surv(RFS_MONTHS,
                    RFS_STATUS) ~ Risk_Ma, km_data_all)
summary(cph_Ma)

fit_Ma <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ Risk_Ma,
                  km_data_all)

sd_Ma <- survdiff(Surv(RFS_MONTHS,
                       RFS_STATUS ) ~ Risk_Ma, km_data_all)
p_val_Ma =1 - pchisq(sd_Ma$chisq, length(sd_Ma$n) - 1)

ggsurvplot(
  fit = fit_Ma, 
  palette = c("red", "black"),
  xlab = "Months", ylab = "Relapse Free Survival",
  title = "Metabric ER+ HT treated patients", subtitle="Ma et al - HOXB13:IL17RB ratio",  
  risk.table=T, tables.y.text=F, 
  pval=p_val_Ma, break.time.by=24,
  ggtheme=theme_survminer(font.main = 20, font.x = 20, font.y = 20, 
                          font.submain = 15, font.tickslab = 20, 
                          font.legend = 20))

DXY_Ma<- cph(Surv(RFS_MONTHS, RFS_STATUS) ~ Risk_Ma,
             km_data_all, x=T, y=T)
validate(DXY_Ma, method='boot', B=150)

##### Cox PH p_val comparison
Cox_ph <- data.frame(Method = factor(c( "5 Cand\nPath", "10GS-Men","HOXB13/\nIL17BR",
                                        "CRISPR","SET\nER/PR", "ODX", "6-GS"),
                                     levels = c( "5 Cand\nPath", "10GS-Men","HOXB13/\nIL17BR",
                                                 "CRISPR", "SET\nER/PR", "ODX", "6-GS")),
                     p_value = c(p_val_Raheem, p_val_Men,p_val_Ma, p_val_Crispr,
                                 p_val_SET, p_val_ODX,p_val_Ours),
                     label = c(0.57, 0.021, 0.012, 1.9*10^-4,1.1*10^-4,9.4*10^-5,
                               6.7*10^-5))

ggplot(Cox_ph, aes(x=Method, y=log10(p_value))) + 
  geom_segment( aes(x=Method, xend=Method, y=-5, yend=log10(p_value), color=Method),
                size=5) +
  geom_label(aes(x=Method, y=log10(p_value)+0.4, label = label,fill=Method)) +
  ggtitle("Cox Proportional Hazard p-values in METABRIC ER+ cohort") +
  theme_linedraw() + guides(color="none", fill="none")+
  theme(panel.grid = element_blank(), text = element_text(size = 20),
        panel.background = element_blank())


