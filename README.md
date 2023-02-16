# BLR-SA-TamResistance


_Filtering.R_ performs DEA analysis and generates the joint set of genes common to the  sequencing datasets read from data. In the example, TCGA breast cancer Tam-treated and MCF7 TamResistant cells. _BLR-SA_routine.R_ implements the SA-like algorithm that tests potential gene signatures from the set of common genes by fitting a Stan GLM model and performing a Metropolis test to keep the most promising genes. It is kept as an independent file because for large enough sets of potential candidates genes it is advisible to run it separatelyon an HPC equipment. From its output, _SA_output_analysis.R_ extracts the final signature. _Validation.R_ shows Survival Analysis plots and fits a Cox Proportional Hazards Model for the selected signature and other major signatures in HT resistance for comparison

For the time-line classification step, divide TCGA's clinical information into separete files for each clinical subcategory (dates, follow-up, treatment info...) as in the example and run _Patient_Classification_Visualization.m_ in Matlab. For an easier visualization in bigger cohorts use _Separating_patients.R_ to generate individual files for good, discards or patients in doubt.



