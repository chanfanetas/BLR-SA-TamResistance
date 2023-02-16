library(rstan)
library(loo)
library(rstanarm)
library(dplyr)


# Functions
Accept_Prob <- function(ACC_current, ACC_new, Temp){
  Prob = min(1, exp(- (ACC_current - ACC_new) / Temp) )
  return(Prob)
}

`%notin%` <- Negate(`%in%`)

#Stan options active
rstan_options(auto_write=T)
options(mc.cores = parallel::detectCores())

#Read data 
data <- read.csv("ALL_HT_log2transformed_Z_score_17genes_wrResponse.csv")
Response <- data[,18]


#Set parameters for SA
n_it = 10000
Temp= 1 # Simulated annealing temperature parameter
k = 2  # SImulated annealing Boltzmann constant equivalent - regulates temperature decrease

#Initialize gene quantities
best_genes <- c(1:17) 
signature <- matrix(nrow = n_it, ncol=27)
elpd_best <- -500
All <- seq(1,17,1)


#SA loop
for (i in 1:n_it){
  
  ## Cooling scheme for SA
  if(i%%2500==0){ 
  Temp=Temp/k
  }
  
  # Selected modification for gene signature
  coin <- runif(1,0,1)
  flag_method = -1 # 1 is subtitution, 2 is removal, 3 is addition
  
  if (length(best_genes)==17){
    flag_method= 2 # Removal
  }else if(length(best_genes)==1){
    flag_method= 3 # Addition
  }else{
    if(coin < 1/3){
      flag_method = 1 # Substitution
    }else if(coin > 2/3){
      flag_method = 2 # Removal
    }else{
      flag_method = 3 # Addition
    }
  }
  
  if(flag_method==1){ 
    # Substitute a gene and create gene score from Z-scores
    rand_from_sig <-  sample(best_genes, length(best_genes)-1)
    All_but_sig <- All[All %notin% best_genes]
    rand_genes <- append(rand_from_sig, sample(All_but_sig,1))
    data_rand <- data.frame(Mean_Rand = apply(data[,rand_genes],1,mean))
    data_rand$Response <- Response
    colnames(data_rand)[1] <- "Gene_Score"
    
    # Fit BLR model
    fit <- stan_glm(formula = Response ~ Gene_Score,
                     data = data_rand, family = binomial(link = "logit"),
                     refresh = 0)
    loo_fit <- loo(fit, save_psis = T)
    preds <- posterior_predict(fit)
    # LOO predictive probabilities
    ploo=E_loo(preds, loo_fit$psis_object, type="mean", log_ratios = -log_lik(fit))$value
    
    # LOO classification accuracy
    TP <- 0 # True positives (are resistant, predicted to be resistant)
    TN <- 0 # True Negatives (Good responders, predicted as GR)
    FP <- 0 # False Positives (are resistant, predicted as GR)
    FN <- 0 # False Negatives (are good responders, predicted to be resistant)
    Sensitivity <- 0
    Specificity <- 0
    Overall <- 0
    for(j in 1:length(ploo)){
      if(ploo[j] < 0.5){ # If we predict GR
        if(Response[j]==0){
          TN <- TN + 1  # Well predicted GR
        }else{
          FP <- FP + 1 # Poorly predicted GR
        }
      }else{ # If we predict Res
        if(Response[j]==1){ 
          TP <- TP + 1 # Well predicted Res
        }else{
          FN <- FN + 1 # Poorly predicted Res
        }
        
      }
    }
    
    Sensitivity <- TP/(TP+FN)
    Specificity <- TN/(TN+FP)
    Overall <- (TP+TN)/(TP+FP+FN+TN)
    
    signature[i,1:length(rand_genes)] <- rand_genes
    signature[i,18] <- loo_fit[[1]][1]
    signature[i,19] <- Sensitivity
    signature[i,20] <- Specificity
    signature[i,21] <- Overall
    signature[i,22] <- TP
    signature[i,23] <- TN
    signature[i,24] <- FP
    signature[i,25] <- FN
    signature[i,26] <- elpd_best
    
    Acceptance <- Accept_Prob(elpd_best, loo_fit[[1]][1], Temp)
    Metropolis_randnumber <- runif(1,0,1)
    if(Acceptance > Metropolis_randnumber){
      best_genes <- rand_genes
      elpd_best <- loo_fit[[1]][1]
      signature[i,27] <- Acceptance
    } 
    
  }else if(flag_method==2){
    # Remove a gene and create gene score from Z-scores
    gene_to_remove <-  sample(length(best_genes), 1) #select an index from the sig length to remove
    rand_genes <- best_genes[-gene_to_remove]
    
    #Fix to calculate mean value when the gene list if 1 gene long
    if(length(rand_genes)==1){
      rand_genes <-rep(rand_genes,2)
    }
    
    
    data_rand <- data.frame(Mean_Rand = apply(data[,rand_genes],1,mean))
    data_rand$Response <- Response
    colnames(data_rand)[1] <- "Gene_Score"
    
    # Fit BLR model
    fit <- stan_glm(formula = Response ~ Gene_Score,
                     data = data_rand, family = binomial(link = "logit"),
                     refresh = 0)
    loo_fit <- loo(fit, save_psis = T)
    preds <- posterior_predict(fit)
    # LOO predictive probabilities
    ploo=E_loo(preds, loo_fit$psis_object, type="mean", log_ratios = -log_lik(fit))$value
    
    # LOO classification accuracy
    TP <- 0 # True positives (are resistant, predicted to be resistant)
    TN <- 0 # True Negatives (Good responders, predicted as GR)
    FP <- 0 # False Positives (are resistant, predicted as GR)
    FN <- 0 # False Negatives (are good responders, predicted to be resistant)
    Sensitivity <- 0
    Specificity <- 0
    Overall <- 0
    for(j in 1:length(ploo)){
      if(ploo[j] < 0.5){ # If we predict GR
        if(Response[j]==0){
          TN <- TN + 1  # Well predicted GR
        }else{
          FP <- FP + 1 # Poorly predicted GR
        }
      }else{ # If we predict Res
        if(Response[j]==1){ 
          TP <- TP + 1 # Well predicted Res
        }else{
          FN <- FN + 1 # Poorly predicted Res
        }
        
      }
    }
    
    Sensitivity <- TP/(TP+FN)
    Specificity <- TN/(TN+FP)
    Overall <- (TP+TN)/(TP+FP+FN+TN)
    
    signature[i,1:length(rand_genes)] <- rand_genes
    signature[i,18] <- loo_fit[[1]][1]
    signature[i,19] <- Sensitivity
    signature[i,20] <- Specificity
    signature[i,21] <- Overall
    signature[i,22] <- TP
    signature[i,23] <- TN
    signature[i,24] <- FP
    signature[i,25] <- FN
    signature[i,26] <- elpd_best

    
    Acceptance <- Accept_Prob(elpd_best, loo_fit[[1]][1], Temp)
    Metropolis_randnumber <- runif(1,0,1)
    if(Acceptance > Metropolis_randnumber){
      best_genes <- rand_genes
      elpd_best <- loo_fit[[1]][1]
      signature[i,27] <- Acceptance
    }


  }else if(flag_method==3){ 
    # Add a gene and create gene score from Z-scores
    All <- seq(1,17,1)
    All_but_sig <- All[All %notin% best_genes]
    rand_genes <- append(best_genes, sample(All_but_sig,1))
    data_rand <- data.frame(Mean_Rand = apply(data[,rand_genes],1,mean))
    data_rand$Response <- Response
    colnames(data_rand)[1] <- "Gene_Score"
    
    # Fit BLR model
    fit <- stan_glm(formula = Response ~ Gene_Score,
                     data = data_rand, family = binomial(link = "logit"),
                     refresh = 0)
    loo_fit <- loo(fit, save_psis = T)
    preds <- posterior_predict(fit)
    # LOO predictive probabilities
    ploo=E_loo(preds, loo_fit$psis_object, type="mean", log_ratios = -log_lik(fit))$value
    
    # LOO classification accuracy
    TP <- 0 # True positives (are resistant, predicted to be resistant)
    TN <- 0 # True Negatives (Good responders, predicted as GR)
    FP <- 0 # False Positives (are resistant, predicted as GR)
    FN <- 0 # False Negatives (are good responders, predicted to be resistant)
    Sensitivity <- 0
    Specificity <- 0
    Overall <- 0
    for(j in 1:length(ploo)){
      if(ploo[j] < 0.5){ # If we predict GR
        if(Response[j]==0){
          TN <- TN + 1  # Well predicted GR
        }else{
          FP <- FP + 1 # Poorly predicted GR
        }
      }else{ # If we predict Res
        if(Response[j]==1){ 
          TP <- TP + 1 # Well predicted Res
        }else{
          FN <- FN + 1 # Poorly predicted Res
        }
        
      }
    }
    
    Sensitivity <- TP/(TP+FN)
    Specificity <- TN/(TN+FP)
    Overall <- (TP+TN)/(TP+FP+FN+TN)
    
    signature[i,1:length(rand_genes)] <- rand_genes
    signature[i,18] <- loo_fit[[1]][1]
    signature[i,19] <- Sensitivity
    signature[i,20] <- Specificity
    signature[i,21] <- Overall
    signature[i,22] <- TP
    signature[i,23] <- TN
    signature[i,24] <- FP
    signature[i,25] <- FN
    signature[i,26] <- elpd_best


    Acceptance <- Accept_Prob(elpd_best, loo_fit[[1]][1], Temp)
    Metropolis_randnumber <- runif(1,0,1)
    if(Acceptance > Metropolis_randnumber){
      best_genes <- rand_genes
      elpd_best <- loo_fit[[1]][1]
      signature[i,27] <- Acceptance
    }


  }
}
# Rename output columns
colnames(signature) <- c(colnames(data[1:17]), "LOO_ELPD", "Sensitivity", "Specificity", "Accuracy","TP","TN","FP","FN","Current_ELPD", "Acceptance")
#  Write result to file
write.csv(signature, "../../data/SA_Halved_Cooling_nit_10K_GSS.csv")
