# ADD-SIS for Interval-Valued Response

## Install Packages: 

`if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")`
    
`BiocManager::install("Icens")`

`install.packages("interval")`

`install.packages("energy")`

## Case Study Folder: 

### Folder Code Description:


| File Name | Description  |
| ------- | --- |
| CaseStudyCode_Final.Rmd | The Rmd file for the case study|

### Folder Data Description:

| File Name | Description  |
| ------- | --- |
| REAL_DATA.csv | The original text data |
|Words.csv| The 1043 words generated from the original text data|
|Words_Filter.csv| The 677 words generated from the original text data by filtering out the words whose frequencies are less than 5|
|WordMatrix_Filter.csv|The big design matrix (X) table consists of 0,1|
|SelectedWord | The several top words based on ADD-SIS associated with the English translation|
______________________________________________________________________________________________
Note: since the data contain Chinese characters, so we need to encode the CSV file as UTF-8. To display these CSV files in the correct way, go to Data > From Text to launch a Text Import Wizard. Now select the file origin to pick ¡°65001: Unicode (UTF-8)¡±, this will turn the CSV files into something that is legible, then choose New Sheet to display the CSV files in the correct way.



## Simulation Folder

### Example1-3  Description:
| File Name | Description  |
| ------- | --- |
| Code/COX_FullVersion.R | The original code to run the simulation under Cox Model setup |
|Code/GLM_FullVersion.R| The original code to run the simulation under Generalized Linear Model  setup |
|Code/LM_FullVersion.R| The original code to run the simulation under  Linear Model  setup|
|Code/ResultsAnalysis.Rmd| Analyze the intermediate results and generate Table 1 in paper|
|Folder Intermediate|Store all the intermediate data for simulations|



### Example 4  Description:
| File Name | Description  |
| ------- | --- |
| Code/Robustness.R | The original code to run the simulation |
|Code/ResultsAnalysis.Rmd| Analyze the intermediate results and generate Table 2 in paper |
|Folder Intermediate| Store all the intermediate data for simulations. |
|Code/ResultsAnalysis.Rmd| Analyze the intermediate results and generate Table 1 in paper|
|Folder Intermediate|Store all the intermediate data for simulations|

# Instructions: 

The following steps give the instruction on how to reproduce the analytical results: 

### Step 1: Install the corresponding packages:

### Step 2: Load the following R- packages with the specified versions: 
| Package | Version  |
| ------- | --- |
|Icens version | 1.58.0|
|Interval version|1.1-0.1|
|MASS version|7.3-51.6|
|dplyr version|0.8.5|
|energy version|1.7-7|
|readr version|1.3.1|


### Step 3: Simulation Example1-3 (Table 1 in paper): Use the “”Simulation/Example1-3” folder:
(1)	Run the “Code/COX_FullVersion.R”, “Code/GLM_ FullVersion.R” and “Code/LM_ FullVersion.R” code.   
(2)	We store the intermediate results for simulation in “Intermediate” folder.  
(3)	Use ‘Code/ResultsAnalysis.Rmd’ to generate performance evaluation for simulation.  

### Step 4: Simulation Exampe 4 (Table 2 in paper):: Use the “”Simulation/Example4” folder:
(1)	Run the “Code/Robustness.R”. It might take 10 hours to run this file.  
(2)	We store the intermediate results for simulation in “Intermediate” folder.  
(3)	Use ‘Code/ResultsAnalysis.Rmd to generate performance evaluation for simulation.  

### Step 5: Case Study Part: Use the “”CaseStudy” folder:
(1)	Use ‘CaseStudy/Code/CaseStudyCode_Final.Rmd’ file to generate the results.  

### Step 6: Supplementary Part: Use the ‘”SUPP” folder: 
(1)	‘TableS.1/TableS.1.Rmd’ reproduces the results in Table S.1.   
(2)	For reproducing the results TableS.2, file ‘TableS.2/Code/DistanceComparison.Rmd’ runs the simulation examples and we store the intermediate results in ‘TableS.2/Intermediate’ folder. File ‘TableS.2/Code/ResultsAnalysis.Rmd’ analyzes the intermediate results.     
(3)	For Table S.3, ‘TableS.3/Overlap.Rmd’ reproduces the results.    
(4)	For Table S.4 and TableS.5, the utilized data is in ‘TableS.4&5/Data’ folder, and ‘TableS.4&5/Code/CaseStudy-Web.Rmd’ runs the corrosponding analysis and the generated results are attahced in ‘TableS.4&5/Results/ CaseStudy-Web.html’ file.   
(5)	For Table S.6-S.8, use the code within ‘TableS.6-8’ folder.  
