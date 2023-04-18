# BLR-SA-TamResistance

Source code for the use and replication of the methodology described in the paper **A novel mathematical approach for analysis of integrated cell-patient data uncovers a 6-gene signature linked to endocrine therapy resistance **. 

- Data

High-throughput data for comparable sets of cells and patients in terms of the biological problem in question is required. Besides, clinical data for said patients, especially concerning follow-up, treatment and response data should be used to classify patients for DEA. The rest of the data files can be derived from them or obtained as the output of some part of the methodology. The examples provided come from already cleaned data from the TCGA and in-lab RNA-seq experiments. Full raw data and GSEA provided in the paper can be found at https://vivancoslabseq.shinyapps.io/RNASeqSOX/.

- Code

_Filtering.R_ performs DEA analysis and generates the joint set of genes common to the sequencing datasets read from the data folder. In the example, TCGA breast cancer Tam-treated and MCF7 TamResistant cells were used. _BLR-SA_routine.R_  implements the SA-like algorithm that tests potential gene signatures from the set of common genes by fitting a Stan GLM model and performing a Metropolis test to keep the most promising genes. It is kept as an independent file because for large enough sets of potential candidate genes it is advisable to run it separately on HPC equipment. From its output, _SA_output_analysis.R_ extracts the final signature. _Validation.R_ Analysis plots and fits a Cox Proportional Hazards Model for the selected signature and other major signatures in HT resistance for comparison.

For the timeline classification step, divide TCGA's clinical information into separate files for each clinical subcategory (dates, follow-up, treatment info...) as in the example and run _Patient_Classification_Visualization.m_ in Matlab. For easier visualization in bigger cohorts use _Separating_patients.R_ to generate individual files for good, discards or patients in doubt.



