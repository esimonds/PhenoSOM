# PhenoSOM


<img src="https://raw.githubusercontent.com/esimonds/PhenoSOM/main/FR-FCM-Z3HK_demo/FR-FCM-Z3HK_Demo_Step2_success.png" alt="tSNE plot of 3900 SOM nodes colored by PhenoGraph cluster" width="250"/>


## What is PhenoSOM?
PhenoSOM is an R script to bring together the published FlowSOM, PhenoGraph, and edgeR algorithms into a single pipeline for differential analysis of cluster abundance in mass cytometry data.

It was used to perform the clustering, differential analyses, and generate the volcano plots in [Simonds et al. _J Immunother Cancer_ (2021)](http://doi.org/10.1136/jitc-2020-002181).


## Should I use PhenoSOM to analyze my own mass cytometry data?
Probably not, for the reasons below. I would recommend using more recent or more actively developed mass cytometry data analysis tools. Check out [Cytoforum](http://cytoforum.stanford.edu) for an active community of mass cytometry users and discussions on the latest algorithms. This script was born in 2016 and many new analysis tools have been published in the years since.

If you're savvy with R, then sure, have fun with it. Just be aware that PhenoSOM doesn't do anything particularly novel -- it simply strings together a few existing R packages that were created by and are maintained by others. PhenoSOM is not an R package. It's just a script, and it is not actively maintained or supported. The help documentation is essentially this README plus the comments in the script itself. I am mainly depositing it here on GitHub for posterity. 


## Who should use PhenoSOM?
The main audience for this GitHub project is bioinformatics researchers who would like to compare and contrast approaches for clustering and/or differential analysis. I think the core workflow of PhenoSOM (e.g. FlowSOM of individual files followed by PhenoGraph of the aggregate) has merit and if it stands up to some more rigorous testing, perhaps someone would be motivated to implement it in a smarter, more user-friendly way. This script evolved a lot over 5 years, and I am not a professional programmer, so I apologize in advance that the code is inelegant, to say the least.

In addiiton to bioinformatics researchers, there is a small subset of mass cytometry users at UCSF that have been using the PhenoSOM script in their own research since 2016, and this project is also meant to serve as a resource for them.



## I'm a bioinformatics researcher and I want to reproduce Figure 1B from Simonds et al (2021) -- how can I do that?
Great! There is a dedicated helper script to make this easier. Please check out the [FR-FCM-Z3HK_demo README](FR-FCM-Z3HK_demo/DemoReadme.md)


## What are the general steps to configure and run PhenoSOM?
1. Install Rstudio (this is required)
2. Install the required packages from BioConductor and CRAN  
```R
install.packages(c("Rtsne", "ggplot2", "RColorBrewer", "gplots", "tidyr", "data.table"))
```
3. Install the required packages from BioConductor and CRAN  
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("flowCore")
BiocManager::install("FlowSOM")

```
4. Install the cytofkit2 package from GitHub **-- not BioConductor! --** as described on the [CytofKit2 GitHub page](https://github.com/JinmiaoChenLab/cytofkit2)
5. Download the R scripts, TXT files, and CSV files from this GitHub project to your analysis folder.
6. Create a subfolder under your analysis folder called **FCSfiles** and put your input files there. It's a good idea to gate the files first to enrich for viable singlets, but the script has been used successfully on ungated files in the past. All input files must have the same staining panel (e.g. set of channels and markers)
7. Modify **clustering markers.txt** to include the channel and marker names of the markers you want to use for clustering. Omit any very noisy markers or general channels (e.g. DNA). Try to include only markers that you think will help identify consistent cell phenotypes across your input files. Be careful including dynamic markers that are not tied to cell fate (e.g. Ki67) because the algorithm is likely to bisect a cell type (e.g. CD8 T cells) based on Ki67-high and Ki67-low (however, in some cases, that's what you want).
8. Modify **plotting markers.txt** to include the channel and marker names of the markers you want to use for plotting. Anything on this list will be included in the generated heatmaps. If you excluded Ki-67 in the step above, you can include it here and it will not affect clustering, but you can see it on the heatmaps.
9. Modify **edgeR_setup.csv** to define a "Condition" associated with each of your input files. The goal here is to do an A vs. B comparison, so you only want two conditions. If you have files that are not part of condition A or B, assign the condition "NA". Note that the Condition value that comes earliest in the alphabet will be used as the denominator or baseline. For example, if you assign the Conditions "DrugA" and "DrugB", the script will express the data as fold-change in the abundance of populations in your "DrugB" samples relative to your "DrugA" samples.
10. Open the **PhenoSOM Step 1.R** script in Rstudio and update the variables under "User-defined parameters" section. Then, run it.
11. Open the **PhenoSOM Step 2.R** script in Rstudio and update the variables under "User-defined parameters" section. Then, run it.
12. Open the **PhenoSOM Step 3.R** script in Rstudio and update the variables under "User-defined parameters" section. Then, run it.


## Something isn't working right. What should I try to troubleshoot it?
Try following the [FR-FCM-Z3HK_demo](FR-FCM-Z3HK_demo/DemoReadme.md) to help determine if it's something with your setup (e.g. packages) or the input files. Once you get the demo working properly and you have a sense of how the different scripts relate to each other, it will be easier to replace the demo files with your own data.


## What package versions are required?
The demo script was tested on a Mac with Rstudio version 1.3.959, R version 4.0.3 and the following package versions:
```R
> packageVersion("Rtsne")
[1] ‘0.15’
> packageVersion("ggplot2")
[1] ‘3.3.3’
> packageVersion("RColorBrewer")
[1] ‘1.1.2’
> packageVersion("gplots")
[1] ‘3.1.1’
> packageVersion("tidyr")
[1] ‘1.1.3’
> packageVersion("data.table")
[1] ‘1.14.0’
> packageVersion("flowCore")
[1] ‘2.0.1’
> packageVersion("FlowSOM")
[1] ‘1.20.0’
> packageVersion("cytofkit2")
[1] ‘2.0.1’
```

## My package versions look fine and I'm still having trouble. Who can I contact for help?
Please file a bug report under https://github.com/esimonds/PhenoSOM/issues and I will get an alert. I will try to respond to new issues with 7 days. I no longer work at UCSF where this script was developed, so maintenance/troubleshooting is something that I have to do in my free time. That being said, I would like others to explore the core worfklow of PhenoSOM and hopefully someone will find it useful, or it will inspire a new, better approach for analyzing mass cytometry data.
