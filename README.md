# PhenoSOM: An R workflow for clustering and differential abundance calculation in mass cytometry data


<img src="https://raw.githubusercontent.com/esimonds/PhenoSOM/main/FR-FCM-Z3HK_demo/images/FR-FCM-Z3HK_Demo_Step2_output_success.png" alt="tSNE plot of 3900 SOM nodes colored by PhenoGraph cluster" width="250"/>


## Q: What is PhenoSOM?
PhenoSOM is a set of R scripts that link together the published FlowSOM, PhenoGraph, and edgeR algorithms into a single pipeline for differential analysis of cluster abundance in mass cytometry data. It is designed to perform "A versus B" comparisons between conditions (e.g. drug vs. vehicle) with multiple replicate samples per condition. It has been used successfully on immunology-focused datasets (e.g. blood, splenocytes, lymph nodes, dissociated tumors). It works well on datasets containing about 10 million cells and about 40 FCS files; it has not been tested on much larger datasets.

The PhenoSOM workflow was used to perform the clustering, differential analyses, and generate the volcano plots in [Simonds et al. _J Immunother Cancer_ (2021)](http://doi.org/10.1136/jitc-2020-002181).


## Q: Should I use PhenoSOM to analyze my own mass cytometry data?
Probably not, for three reasons. Firstly, this script was born in 2016 and many new analysis tools have been published in the years since. I would recommend using more recent or more actively developed mass cytometry data analysis tools. Check out [Cytoforum](http://cytoforum.stanford.edu) for an active community of mass cytometry users and discussions on the latest algorithms.

Secondly, PhenoSOM doesn't do anything particularly novel -- it simply strings together a few existing R packages that were created by and are maintained by others. PhenoSOM is not an R package. It's just a set of scripts, and it is not actively maintained or supported. 

Thirdly, the help documentation is limited to what you find in this README, [the Wiki](https://github.com/esimonds/PhenoSOM/wiki), and in the comments in the script itself. I am mainly depositing it here on GitHub for posterity. 

If you're savvy with R and, despite all these caveats, you want to run your data through PhenoSOM, that's great, and I hope it's useful to you.


## Q: Who should use PhenoSOM?
The main audience for this GitHub project is bioinformatics researchers who would like to compare and contrast approaches for clustering and/or differential analysis. I think the core workflow of PhenoSOM is powerful and useful (e.g. FlowSOM of individual files followed by PhenoGraph of the aggregate; then repeating this on the CD45-positive cells). Perhaps, if it stands up to some more rigorous testing, someone would be motivated to implement it in a smarter, more user-friendly way. This script evolved a lot over 5 years, and I am not a professional programmer, so I apologize that the code is inelegant, to say the least.

In addition to bioinformatics researchers, there is a small subset of mass cytometry users at UCSF that have been using the PhenoSOM script in their own research since 2016, and this project is also meant to serve as a resource for them.



## Q: How can I get started running PhenoSOM?
The best way to get started with PhenoSOM is to run it on the demo dataset. There is a dedicated helper script to make this easier. First, follow the installation instructions below, and then visit [Running the FR-FCM-Z3HK_demo](https://github.com/esimonds/PhenoSOM/wiki/Running-the-FR-FCM-Z3HK-demo)


## Q: How do I install PhenoSOM?
PhenoSOM is a set of R scripts -- it is not an R package that requires installation. You need to copy the scripts to a folder containing your data (FCS files) and the required configuration files (CSV and TXT files). However, several R packages are required for PhenoSOM to work. Visit the [Installing PhenoSOM](https://github.com/esimonds/PhenoSOM/wiki/Installing-PhenoSOM) for details on how to install those.


## Q. How do I configure and run PhenoSOM on my own data?
First, install the required R packages as described above. Perferably, follow the [FR-FCM-Z3HK_demo](https://github.com/esimonds/PhenoSOM/wiki/Running-the-FR-FCM-Z3HK-demo) to make sure everything is working properly and learn how PhenoSOM's configuration files and scripts work together. Then, visit [Configuring PhenoSOM](https://github.com/esimonds/PhenoSOM/wiki/Configuring-PhenoSOM) for details on how to configure the scripts for your own dataset.


## Q: Something isn't working right. What should I try to troubleshoot it?
Try following the [FR-FCM-Z3HK_demo](https://github.com/esimonds/PhenoSOM/wiki/Running-the-FR-FCM-Z3HK-demo) to help determine if it's something with your system (e.g. R packages) or with your PhenoSOM configuration. Once you get the demo working properly and you have a sense of how the different scripts relate to each other, it will be easier to replace the demo files with your own data.


## Q: What package versions are required?
More info at [Installing PhenoSOM](https://github.com/esimonds/PhenoSOM/wiki/Installing-PhenoSOM)


## Q: My package versions look fine and I'm still having trouble. Who can I contact for help?
Please file a bug report under https://github.com/esimonds/PhenoSOM/issues and I will get an alert. I will try to respond to new issues with 7 days. I no longer work at UCSF where this script was developed, so maintenance/troubleshooting is something that I have to do in my free time. That being said, I would like others to explore the core worfklow of PhenoSOM and hopefully someone will find it useful, or it will inspire better approaches for analyzing mass cytometry data.
