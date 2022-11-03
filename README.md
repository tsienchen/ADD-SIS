# ADD-SIS

## Install Packages: 

`if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")`
`BiocManager::install("Icens")`

`install.packages(“Interval”)`


## Case Study 

### Folder Code Description:

CaseStudyCode_Final.R: The R code for the case study. 



### Folder Data Description:

REAL_DATA.csv: The original text data 

Words.csv: 1043 words generated from the original text data

Words_Filter.csv: 677 words generated from the original text data by filtering out the words whose frequencies are less than 5.

WordMatrix_Filter.csv: The big design matrix (X) table consists of 0,1

SelectedWord: The several top words based on ADD-SIS associated with the English translation


______________________________________________________________________________________________
Note: since the data contain Chinese characters, so we need to encode the CSV file as UTF-8. To display these CSV files in the correct way, go to Data > From Text to launch a Text Import Wizard. Now select the file origin to pick ¡°65001: Unicode (UTF-8)¡±, this will turn the CSV files into something that is legible, then choose New Sheet to display the CSV files in the correct way.



## Simulation

### Example1-3 Folder Description:

Code/COX_FullVersion.R,Code/GLM_FullVersion.R, Code/LM_FullVersion.R: The original code to tun the simulation

Code/ResultsAnalysis.Rmd: Analyze the intermediate results and generate Table 1 in paper. 
Intermediate: Store all the intermediate data for simulations. 
Results/ResultsAnalysis.html: The returned results for Code/ResultsAnalysis.Rmd. 



### Example4 Folder Description:

Code/Robustness.R: The original code to tun the simulation

Code/ResultsAnalysis.Rmd: Analyze the intermediate results and generate Table 2 in paper. 
Intermediate: Store all the intermediate data for simulations. 
Results/ResultsAnalysis.html: The returned results for Code/ResultsAnalysis.Rmd. 
