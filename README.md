# PhenoSOM: An R workflow for clustering and differential abundance calculation in mass cytometry data


<img src="https://raw.githubusercontent.com/esimonds/PhenoSOM/main/FR-FCM-Z3HK_demo/images/FR-FCM-Z3HK_Demo_Step2_output_success.png" alt="tSNE plot of 3900 SOM nodes colored by PhenoGraph cluster" width="250"/>



# Note:  Code is not yet visible in this repository. I'm just getting the documentation up first. The code will be posted when the article is published.



## Q: What is PhenoSOM?
PhenoSOM is an R script that links the published FlowSOM, PhenoGraph, and edgeR algorithms into a single pipeline for differential analysis of cluster abundance in mass cytometry data.

It was used to perform the clustering, differential analyses, and generate the volcano plots in [Simonds et al. _J Immunother Cancer_ (2021)](http://doi.org/10.1136/jitc-2020-002181).


## Q: Should I use PhenoSOM to analyze my own mass cytometry data?
Probably not, for the reasons below. I would recommend using more recent or more actively developed mass cytometry data analysis tools. Check out [Cytoforum](http://cytoforum.stanford.edu) for an active community of mass cytometry users and discussions on the latest algorithms. This script was born in 2016 and many new analysis tools have been published in the years since.

Please also be aware that PhenoSOM doesn't do anything particularly novel -- it simply strings together a few existing R packages that were created by and are maintained by others. PhenoSOM is not an R package. It's just a script, and it is not actively maintained or supported. The help documentation is essentially this README plus the comments in the script itself. I am mainly depositing it here on GitHub for posterity. If you're savvy with R and, despite all these caveats, you want to run your data through PhenoSOM, then have fun and I hope it's useful to you.


## Q: Who should use PhenoSOM?
The main audience for this GitHub project is bioinformatics researchers who would like to compare and contrast approaches for clustering and/or differential analysis. I think the core workflow of PhenoSOM (e.g. FlowSOM of individual files followed by PhenoGraph of the aggregate) has merit and if it stands up to some more rigorous testing, perhaps someone would be motivated to implement it in a smarter, more user-friendly way. This script evolved a lot over 5 years, and I am not a professional programmer, so I apologize in advance that the code is inelegant, to say the least.

In addiiton to bioinformatics researchers, there is a small subset of mass cytometry users at UCSF that have been using the PhenoSOM script in their own research since 2016, and this project is also meant to serve as a resource for them.



## Q: How can I get started running PhenoSOM?
The best way to get started with PhenoSOM is to run it on the demo dataset. There is a dedicated helper script to make this easier. First, follow the installation instructions below, and then visit [Running the FR-FCM-Z3HK_demo](https://github.com/esimonds/PhenoSOM/wiki/Running-the-FR-FCM-Z3HK-demo)


## Q: How do I install PhenoSOM?
PhenoSOM is a series of R scripts -- it is not an R package. However, some existing R packages are required for PhenoSOM to work. Visit the [Installing PhenoSOM](https://github.com/esimonds/PhenoSOM/wiki/Installing-PhenoSOM) for details.


## Q. How do I configure and run PhenoSOM on my own data?
First, install the required R packages as described above. Perferably, follow the [FR-FCM-Z3HK_demo](https://github.com/esimonds/PhenoSOM/wiki/Running-the-FR-FCM-Z3HK-demo) to make sure everything is working properly. Then, visit [Configuring PhenoSOM](https://github.com/esimonds/PhenoSOM/wiki/Configuring-PhenoSOM) for details on how to configure the scripts for your own dataset.


## Q: Something isn't working right. What should I try to troubleshoot it?
Try following the [FR-FCM-Z3HK_demo](https://github.com/esimonds/PhenoSOM/wiki/Running-the-FR-FCM-Z3HK-demo) to help determine if it's something with your setup (e.g. packages) or the input files. Once you get the demo working properly and you have a sense of how the different scripts relate to each other, it will be easier to replace the demo files with your own data.


## Q: What package versions are required?
More info at [Installing PhenoSOM](https://github.com/esimonds/PhenoSOM/wiki/Installing-PhenoSOM)


## Q: My package versions look fine and I'm still having trouble. Who can I contact for help?
Please file a bug report under https://github.com/esimonds/PhenoSOM/issues and I will get an alert. I will try to respond to new issues with 7 days. I no longer work at UCSF where this script was developed, so maintenance/troubleshooting is something that I have to do in my free time. That being said, I would like others to explore the core worfklow of PhenoSOM and hopefully someone will find it useful, or it will inspire better approaches for analyzing mass cytometry data.
