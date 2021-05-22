# PhenoSOM FR-FCM-Z3HK Demo
## Dataset:  Simonds et al. 2021 Figure 1A-C; FlowRepository accession ID FR-FCM-Z3HK

Note: The input data for this demo is 2.86GB. It takes about 3 hours (mostly unattended) to run this demo all the way through, from raw data to volcano plots. 

## Goal of the demo:
By following this demo, you will recreate the plots in Figure 1A (center panel) and Figure 1B (volcano plot and heatmap) from Simonds et al *JITC* (2021):

<img src="https://github.com/esimonds/PhenoSOM/raw/main/FR-FCM-Z3HK_demo/images/Simonds_et_al_JITC_2021_Figure_1.png" alt="Simonds et al JITC 2021 Figure 1" width="250"/>

## Steps to run the demo:

1. Install the required PhenoSOM packages as described on the [main readme](../README.md)
2. Download and unzip **PhenoSOM demo FR-FCM-Z3HK.zip** from [here](https://github.com/esimonds/PhenoSOM/raw/main/FR-FCM-Z3HK_demo/PhenoSOM_demo_FR-FCM-Z3HK.zip)
3. Visit [FlowRepository accession ID FR-FCM-Z3HK](http://flowrepository.org/experiments/3636) and download the FCS files from the experiment with "Tcellpanel" in the filename. There should be 39 FCS files totaling about 2.86GB. Move these to the **FCSfiles** subfolder.
4. Run the **PhenoSOM demo FR-FCM-Z3HK master script.R** script and follow the prompts
5. After running Step 1, check if your output matches the figure in the demo. Open the PNG file located at **Your_analysis_folder/Maps and plots of SOM nodes/FR-FCM-Z3HK Panel 1 PhenoSOM Step 1 tSNE map of SOM nodes by Rphenograph metacluster k_30 and size scaled to cell number.png** 

It should look like this:

<img src="https://github.com/esimonds/PhenoSOM/raw/main/FR-FCM-Z3HK_demo/images/FR-FCM-Z3HK_Demo_Step1_output_success.png" alt="tSNE plot of Step 1 SOM nodes colored by PhenoGraph cluster" width="250"/>

6. After running Step 2, check if your output matches the figure in the demo. Open the PNG file located at **Your_analysis_folder/Maps and plots of SOM nodes/FR-FCM-Z3HK Panel 1 PhenoSOM Step 2 tSNE map of SOM nodes by Rphenograph metacluster k_30 and size scaled to cell number.png** 

It should look like this:

<img src="https://github.com/esimonds/PhenoSOM/raw/main/FR-FCM-Z3HK_demo/images/FR-FCM-Z3HK_Demo_Step2_output_success.png" alt="tSNE plot of Step 2 SOM nodes colored by PhenoGraph cluster" width="250"/>

7. Continue with running Step 3.
8. After running Step 3, check if your output matches the figure in the demo. Open the PNG file located at **Your_analysis_folder/PhenoSOM_Step3_output/edgeR_Run1_Glioma_vs_Kidney/edgeR_Run1_Glioma_vs_Kidney metaclusters edgeR volcano plot.pdf**

It should look like this:

<img src="https://github.com/esimonds/PhenoSOM/raw/main/FR-FCM-Z3HK_demo/images/FR-FCM-Z3HK_Demo_Step3_output_success.png" alt="Volcano plot of Glioma vs Kidney" width="250"/>

## Q: How long does each step take?
A:  About 2.5 hours for Step 1, 25 minutes for Step 2, and 1 minute for Step 3 (timing is based on a 2016 Macbook Pro).

## Q:  If the clustering results are stochastic, why does my output look identical to the demo?
A:  The script uses the command "set.seed(42)" to force a random seed of 42. This was the same seed used when running the data in the Simonds et al. (2021) paper and in the demo.
