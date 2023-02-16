# BLR-SA-TamResistance

Source code for the use and repliation of the methodology described in the paper **A novel mathematical approach for the combined analysis of cell and patient data uncovers a 6-gene signature linked to endocrine therapy resistance**. 

- Data
The minimum necessary data is high-throughput data for comparable sets of cell and patients in terms of the biological problem in question, as well as the clinical data for said patients, most importantly, follow-up, treatment and response data. The rest of the data files can be derived from them or obtained as the output of some part of the methodology. The examples provided come from already cleaned data from the TCGA and in lab RNA-seq experients. Full raw data and GSEA provided in the paper can be found at https://vivancoslabseq.shinyapps.io/RNASeqSOX/.

- Code
_Filtering.R_ performs DEA analysis and generates the joint set of genes common to the sequencing datasets read from the _data_ folder. In the example, TCGA breast cancer Tam-treated and MCF7 TamResistant cells were used. _BLR-SA_routine.R_ implements the SA-like algorithm that tests potential gene signatures from the set of common genes by fitting a Stan GLM model and performing a Metropolis test to keep the most promising genes. It is kept as an independent file because for large enough sets of potential candidates genes it is advisible to run it separatelyon an HPC equipment. From its output, _SA_output_analysis.R_ extracts the final signature. _Validation.R_ shows Survival Analysis plots and fits a Cox Proportional Hazards Model for the selected signature and other major signatures in HT resistance for comparison

For the time-line classification step, divide TCGA's clinical information into separete files for each clinical subcategory (dates, follow-up, treatment info...) as in the example and run _Patient_Classification_Visualization.m_ in Matlab. For an easier visualization in bigger cohorts use _Separating_patients.R_ to generate individual files for good, discards or patients in doubt.



