# Drug repositioning for SARS-CoV-2 by Gaussian kernel similarity bilinear matrix factorization .
Coronavirus disease 2019 (COVID-19), a disease caused by severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2), is currently spreading rapidly around the world. Since SARS-CoV-2 seriously threatens human life and health as well as the development of the world economy, it is very urgent to identify effective drugs against this virus. However, traditional methods to develop new drugs are costly and time-consuming, which makes drug repositioning a promising exploration direction for this purpose. In this study, we collected known antiviral drugs to form five virus-drug association datasets, and then explored drug repositioning for SARS-CoV-2 by Gaussian kernel similarity bilinear matrix factorization (VDA-GKSBMF). By the 5-fold cross-validation, we found that VDA-GKSBMF has an area under curve (AUC) value of 0.8851, 0.8594, 0.8807, 0.8824, and 0.8804 respectively on the five datasets, which are higher than those of other state-of-art algorithms in four datasets. Based on known virus-drug association data, we used VDA-GKSBMF to prioritize the top-k candidate antiviral drugs that are most likely to be effective against SARS-CoV-2. We confirmed that the top-10 drugs can be molecularly docked with virus spikes protein/human ACE2 by AutoDock on five datasets. Among them, four antiviral drugs ribavirin, remdesivir, oseltamivir, and zidovudine have been under clinical trials or supported in recent literatures. The results suggest that VDA-GKSBMF is an effective algorithm for identifying potential antiviral drugs against SARS-CoV-2.     


## Requirements
Matlab 2016 or above.   


## Codes 
A_VDA_GMSBMF.m: VDA-GKSBMF algorithm  <br>
% Input:  <br>
matDV: drug-virus matrix <br> 
Wdd: drug-drug matrix <br> 
Wvv: virus-virus matrix <br> 
others: parameters <br> 
% Ouput: <br>
M_recovery: matrix for predicted drug-virus scores <br> 

%main.m: cross-validation code.  <br>


## Dataset
Data is located in the directory: ./VDdataset1 <br>    
matDrugVirus.txt <-- drug-virus matrix   <br> 
matDrugDrug.txt <-- drug-drug matrix <br> 
matVirusVirus.txt <-- virus-virus matrix   <br> 


## Results 
The results will be automatically saved into the directory: Results.   <br>

## cite
If you use this code in your research, please cite: <br> 
Yang, et al. Computational drug repositioning based on multi-similarities bilinear matrix factorization. Briefings in Bioinformatics 22.4 (2021): bbaa267. <br> 
Wang, et al. Drug repositioning for SARS-CoV-2 by Gaussian kernel similarity bilinear matrix factorization. <br>  



## contact<br>
Email: xiang.ju@foxmail.com  
