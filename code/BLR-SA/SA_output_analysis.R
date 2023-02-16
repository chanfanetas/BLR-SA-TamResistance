library(rstan)
library(loo)
library(rstanarm)
library(dplyr)

setwd("~/bayes-cancer/BLR-SA cancer/code/")

#Stan options active
rstan_options(auto_write=T)
options(mc.cores = parallel::detectCores())

#Read data fed to the BLR-SA algorithm and output 
data <- read.csv("../data/ALL_HT_log2transformed_Z_score_17genes_wrResponse.csv")
Response <- data[,18]
data <- data[,1:17]
SA_res<- read.csv("../data/SA_Halved_Cooling_nit_10K_GSS.csv")
SA_res <- SA_res[,-1]

# Remove some vectors with duplicated values (failed suggestions in BLR-SA)
dup <- c()
for(i in 1:dim(SA_res)[1]){
  a <- c()
  a <- as.matrix(SA_res[i,1:17]) %>% as.vector()
  
  if( sum(duplicated(na.omit(a)))>=1){
    dup <- append(dup, i)
  }
}
SA_res_nodup <- SA_res[-dup,]

# Set NA's to zero to avoid errors handling vectors
All_ordered <- matrix(nrow = dim(SA_res)[1], ncol = 17)

for(i in 1:dim(SA_res_nodup)[1]){
  Add_zeros <- c()
  Add_zeros <- as.matrix(SA_res_nodup[i,1:17]) %>% as.vector()
  Add_zeros <- Add_zeros[order(Add_zeros)]
  Add_zeros[is.na(Add_zeros)] <- 0
  All_ordered[i,] <- Add_zeros
}

to_remove <- c()
flag_next <- 0
for(i in 2:dim(SA_res_nodup)[1]){
  for (j in 1:(i-1)){
    if (flag_next==1){
      flag_next=0
      next()
    }
    if( all(All_ordered[i,] == All_ordered[j,])){
      to_remove <- append(to_remove, i)
      flag_next <- 1
    }
  }
}
to_remove <- unique(to_remove)

SA_res_final <- SA_res_nodup[-to_remove,]

# Select al gene combinations with max Accuracy
Result_ACC <- SA_res_final[which(SA_res_final$Accuracy==max(SA_res_final$Accuracy)),]
# Resolve ties using ELPD
Result_ELPD <- Result_ACC[which(Result_ACC$Current_ELPD==max(Result_ACC$Current_ELPD)),]

# Get gene components of the final signature
Gene_signature <- Result_ELPD[1:17]
Gene_signature <- Gene_signature[!is.na(Gene_signature)]
# Print Gene Name of the genes in the signature
colnames(data)[Gene_signature]

# From the final gene signature, recalculate gene score
data_loo <- data.frame(Mean_Sig = apply(data[,Gene_signature],1,mean))
data_loo$Response <- Response

# Fit BLR model and obtain loo predictions
fit <- stan_glm(formula = Response ~ Mean_Sig,
                 data = data_loo, family = binomial(link = "logit"),
                 refresh = 0)

loo <- loo(fit, save_psis = T)
loo

preds <- posterior_epred(fit)
pred <- colMeans(preds)

# LOO predictive probabilities
ploo=E_loo(preds, loo$psis_object, type="mean",
           log_ratios = -log_lik(fit))$value


######### Plots

# Plot of correlation between ELPD and ACC

SA_res_final$FactorACC <- as.factor(SA_res_final$Accuracy)
SA_res_final$Combination <- rownames(SA_res_final)
ggplot(data=SA_res_final,aes(x=Accuracy, y=LOO_ELPD, group=FactorACC) ) + 
  geom_boxplot(alpha=0.2) + xlab(label = "Prediction Accuracy") +
  ylab(label = "Model ELPD") +  theme_linedraw() +
  theme(panel.grid = element_blank(), text = element_text(size = 20),
        panel.background = element_blank() ) 



### Accuracy plot

Plot_acc <- data.frame(Predicted = ploo, Patients = condition_All$Patient, 
                       True = Response,  color = as.factor(Response))
Bad_responses <- Plot_acc[which(Plot_acc$Predicted<0.15 & Plot_acc$True==1),]
Plot_acc <- Plot_acc[-which(Plot_acc$Patients %in% Bad_responses$Patients),]

ggplot(Plot_acc) + geom_point(aes(x=Predicted, y=True,colour=color), size=2)+
  geom_segment(aes(x = min(ploo), y = min(ploo), xend = max(ploo), yend = max(ploo)), size=1) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limit=c(0, 1)) +
  labs(x="Predicted response", y="True Response") + 
  scale_color_manual(values=c("lightblue", "red"),  name="Patient\nResponse",
                     labels=c("Good", "Resistant")) + theme_linedraw()  +
  theme(panel.grid = element_blank(), text = element_text(size = 20), 
        axis.title =  element_text(size = 22), panel.background = element_blank())


