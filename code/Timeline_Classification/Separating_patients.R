library(tidyverse)

Classification <- readxl::read_excel("Clasifying_TCGA_barcodes.xlsx")

#raw data for extra info
raw <- read.csv(file = "Tamoxifen_treated_patients_clinicalinfo_and_drugschedule.csv")
raw_followup <- read.csv(file = "Tamoxifen_treated_patients_followup_radiation_metastasis.csv")

Dates <- read.csv(file = "ClinicalData_dates.csv")
GenderStatus <- read.csv(file = "ClinicalData_GenderAndStatus.csv")
LastFollowup <- read.csv(file = "ClinicalData_LastFollowup.csv")
Neoplasm <- read.csv(file = "ClinicalData_NeoplasmCancerStatus.csv")
NewTumour <- read.csv(file = "ClinicalData_NewTumor.csv")
Radiation <- read.csv(file = "ClinicalData_RadiationInfo.csv")
TreatmentStatus <- read.csv(file = "ClinicalData_TreatmentStatus.csv")
TreatmentType <- read.csv(file = "ClinicalData_TreatmentType.csv")

#CHECK FOR MISSED ANY BARCODEs (Should return no columns)
Classification[!(Classification$Keep %in% Dates$bcr_patient_barcode),Classification$Keep]
Classification[!(Classification$Discard %in% Dates$bcr_patient_barcode), Classification$Discard]
Classification[!(Classification$`Doubtful treatment` %in% Dates$bcr_patient_barcode),Classification$`Doubtful treatment`]

#### Good Treatment
Good_Dates <- Dates[which(Dates$bcr_patient_barcode %in% Classification$Keep),]
Good_GenderStatus <- GenderStatus[which(GenderStatus$bcr_patient_barcode %in% Classification$Keep),]
Good_LastFollowup <- LastFollowup[which(LastFollowup$bcr_patient_barcode %in% Classification$Keep),]
Good_Neoplasm <- Neoplasm[which(Neoplasm$bcr_patient_barcode %in% Classification$Keep),]
Good_NewTumour <- NewTumour[which(NewTumour$Code %in% Classification$Keep),]
Good_Radiation <- Radiation[which(Radiation$bcr_patient_barcode %in% Classification$Keep),]
Good_TreatmentStatus <- TreatmentStatus[which(TreatmentStatus$bcr_patient_barcode %in% Classification$Keep),]
Good_TreatmentType <- TreatmentType[which(TreatmentType$bcr_patient_barcode %in% Classification$Keep),]

Raw_Good <- raw[which(raw[,1] %in% Classification$Keep),]

write.csv(Good_Dates, file = "Good/GoodClinicalData_dates.csv", row.names = F, na="NaN")
write.csv(Good_GenderStatus, file = "Good/GoodClinicalData_GenderAndStatus.csv",  row.names = F)
write.csv(Good_LastFollowup, file = "Good/GoodClinicalData_LastFollowup.csv",  row.names = F, na="NaN")
write.csv(Good_Neoplasm, file = "Good/GoodClinicalData_NeoplasmCancerStatus.csv",  row.names = F, na="NaN")
write.csv(Good_NewTumour, file = "Good/GoodClinicalData_NewTumor.csv",  row.names = F, na="NaN")
write.csv(Good_Radiation, file = "Good/GoodClinicalData_RadiationInfo.csv",  row.names = F, na= "")
write.csv(Good_TreatmentStatus, file = "Good/GoodClinicalData_TreatmentStatus.csv",  row.names = F)
write.csv(Good_TreatmentType, file = "Good/GoodClinicalData_TreatmentType.csv",  row.names = F)

### Doubtful treatment (needed for plotting visualization and having extra visual check)###
Classification <- Classification[order(Classification$`Why Doubtful`),]

TAM300 <- Classification[1:26,]
TAM500 <- Classification[27:44,]
Chaotic <- Classification[45:54,]
TAM1000 <- Classification[55:60,]
TAMIncomplete <- Classification[61:81,]
ParallelTherapy <- Classification[82:87,]

#Select one or more from above to use for dataframe and group
TAMDoubt <- rbind(TAM500,TAM1000)

Raw_TAM_Doubt <- raw[which(raw[,1] %in% TAMDoubt$`Doubtful treatment`),]
Raw_TAM_Doubt_Followup <- raw_followup[which(raw_followup[,1] %in% TAMDoubt$`Doubtful treatment`),]


Doubt_Dates <- Dates[which(Dates$bcr_patient_barcode %in% TAMDoubt$`Doubtful treatment`),]
Doubt_GenderStatus <- GenderStatus[which(GenderStatus$bcr_patient_barcode %in% TAMDoubt$`Doubtful treatment`),]
Doubt_LastFollowup <- LastFollowup[which(LastFollowup$bcr_patient_barcode %in% TAMDoubt$`Doubtful treatment`),]
Doubt_Neoplasm <- Neoplasm[which(Neoplasm$bcr_patient_barcode %in% TAMDoubt$`Doubtful treatment`),]
Doubt_NewTumour <- NewTumour[which(NewTumour$Code %in% TAMDoubt$`Doubtful treatment`),]
Doubt_Radiation <- Radiation[which(Radiation$bcr_patient_barcode %in% TAMDoubt$`Doubtful treatment`),]
Doubt_TreatmentStatus <- TreatmentStatus[which(TreatmentStatus$bcr_patient_barcode %in% TAMDoubt$`Doubtful treatment`),]
Doubt_TreatmentType <- TreatmentType[which(TreatmentType$bcr_patient_barcode %in% TAMDoubt$`Doubtful treatment`),]

write.csv(Doubt_Dates, file = "Doubt/DoubtClinicalData_dates.csv", row.names = F, na="NaN")
write.csv(Doubt_GenderStatus, file = "Doubt/DoubtClinicalData_GenderAndStatus.csv",  row.names = F)
write.csv(Doubt_LastFollowup, file = "Doubt/DoubtClinicalData_LastFollowup.csv",  row.names = F, na="NaN")
write.csv(Doubt_Neoplasm, file = "Doubt/DoubtClinicalData_NeoplasmCancerStatus.csv",  row.names = F, na="NaN")
write.csv(Doubt_NewTumour, file = "Doubt/DoubtClinicalData_NewTumor.csv",  row.names = F, na="NaN")
write.csv(Doubt_Radiation, file = "Doubt/DoubtClinicalData_RadiationInfo.csv",  row.names = F, na="")
write.csv(Doubt_TreatmentStatus, file = "Doubt/DoubtClinicalData_TreatmentStatus.csv",  row.names = F)
write.csv(Doubt_TreatmentType, file = "Doubt/DoubtClinicalData_TreatmentType.csv",  row.names = F)


#### Discard
Discard_Dates <- Dates[which(Dates$bcr_patient_barcode %in% Classification$Discard),]
Discard_GenderStatus <- GenderStatus[which(GenderStatus$bcr_patient_barcode %in% Classification$Discard),]
Discard_LastFollowup <- LastFollowup[which(LastFollowup$bcr_patient_barcode %in% Classification$Discard),]
Discard_Neoplasm <- Neoplasm[which(Neoplasm$bcr_patient_barcode %in% Classification$Discard),]
Discard_NewTumour <- NewTumour[which(NewTumour$Code %in% Classification$Discard),]
Discard_Radiation <- Radiation[which(Radiation$bcr_patient_barcode %in% Classification$Discard),]
Discard_TreatmentStatus <- TreatmentStatus[which(TreatmentStatus$bcr_patient_barcode %in% Classification$Discard),]
Discard_TreatmentType <- TreatmentType[which(TreatmentType$bcr_patient_barcode %in% Classification$Discard),]

Raw_Discard <- raw[which(raw[,1] %in% Classification$Discard),]

write.csv(Discard_Dates, file = "Discard/DiscardClinicalData_dates.csv", row.names = F, na="NaN")
write.csv(Discard_GenderStatus, file = "Discard/DiscardClinicalData_GenderAndStatus.csv",  row.names = F)
write.csv(Discard_LastFollowup, file = "Discard/DiscardClinicalData_LastFollowup.csv",  row.names = F, na="NaN")
write.csv(Discard_Neoplasm, file = "Discard/DiscardClinicalData_NeoplasmCancerStatus.csv",  row.names = F, na="NaN")
write.csv(Discard_NewTumour, file = "Discard/DiscardClinicalData_NewTumor.csv",  row.names = F, na="NaN")
write.csv(Discard_Radiation, file = "Discard/DiscardClinicalData_RadiationInfo.csv",  row.names = F, na= "")
write.csv(Discard_TreatmentStatus, file = "Discard/DiscardClinicalData_TreatmentStatus.csv",  row.names = F)
write.csv(Discard_TreatmentType, file = "Discard/DiscardClinicalData_TreatmentType.csv",  row.names = F)

