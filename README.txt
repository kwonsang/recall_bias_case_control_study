These files are the Rcode files to apply the estimation methods presented in Lee and Dominici's "Accounting for recall bias in case-control studies: a causal inference approach." 
----------------------------------------------------------------
Overview 

(1) Application
This part is aimed to replicating the results from the Wisconsin Logitudinal Study (WLS) dataset discussed in Section 6. We are unable to share the original data, but it can be downloadable at https://www.ssc.wisc.edu/wlsresearch/data/. We include a Rcode file, called "WLS_to_dataset.R", that can create the same dataset that we used in the paper.  Figures 2, 3 and Tables 3, 4 can be replicated with this dataset. 

(2) Simulation 
This part is aimed to replicating the simulation results such as Figure 1 & Table 2 in the paper.


----------------------------------------------------------------
Description of Directories

There are three directories (i) /Rfunctions, (ii) /Rcode and (iii) /Simulation

--
(i) Rfunctions - contains one R script "basic_functions". This script file contains all the R functions that implement the estimation methods discussed in the paper. The maximum likelihood and stratification methods are implemented. Also, we implement the MH method discussed in the supplementary materials.  

--
(ii) Rcode - contains three R scripts and one csv file. 
	(a) "WLS_to_dataset.R" - this shows how to create a dataset used in Section 6 from the WLS data. 
	(b) "data_analysis.R" - it contains all the statistical analysis. The two estimation methods are implemented here with the dataset obtained from (a). Also, we illustrate how readers can reproduce Tables 3, 4 & Figures 2, 3. 
	(c) "replication.R" - it creates the same tables and figures shown in the manuscript. 
	(d) "sensi_res.csv" - it contains the results of statistical analysis. 

--
(iii) simulation - contains two R scripts and one sub-folder containing simulation results. 
	(a) "simple_recall_bais_numerical_example.R" - shows the simple small simulation study illustrated in Section 3.3. 
	(b) "sim_rep.R" - shows the simulation study discussed in Section 5. 
	(c) "sim_res" folder - contains the simulation results that can be obtained from running sim_rep.R for all different simulation scenarios. 
		(c1) "analysis.R" - contains R codes that can replicate Table 2. 


----------------------------------------------------------------

Computer specification
	Model Name: iMac Pro
	Processor Name: 8-Core Intel Xeon W
	Processor Speed: 3.2 GHz
	Number of Processors:	1
	Total Number of Cores: 8
	Memory: 64 GB
----------------------------------------------------------------
Computing Information:
R version 3.6.3 (2020-02-29)Platform: x86_64-apple-darwin15.6.0 (64-bit)Running under: macOS Catalina 10.15.7


To run the R scripts, several packages should be installed in advance. 

Attached packages:
(Data analysis)
doParallel_1.0.16
iterators_1.0.13
foreach_1.5.1   
optimx_2020-4.2  

(Data mingling)
labelled_2.7.0  
haven_2.3.1  
forcats_0.5.0
stringr_1.4.0
dplyr_1.0.4
purrr_0.3.4    
readr_1.4.0
tidyr_1.1.2
tibble_3.0.4
ggplot2_3.3.2
tidyverse_1.3.0
foreign_0.8-80 
