# Introduction

This document provides information for using the [CoralSNP Galaxy environment](https://coralsnp.science.psu.edu/galaxy)
which is based on the [Galaxy workbench](https://galaxyproject.org/).  A general understanding of Galaxy is required, so please visit
the [introductory tutorials](https://training.galaxyproject.org/training-material/topics/introduction) if you are not yet familiar
with Galaxy.

The Galaxy CoralSNP environment enables streamlined analysis of the coral SNPchip available from Fisher Scientific to ultimately
provide the user with a genet id, converted raw genotyped data, sample relatedness and hybrid status.

The process is straightforward.  Each of the following steps will be discussed in detail in the following sections.

 - A sample metadata file is created by the user from an [Excel spreadsheet template](http://baumslab.org/documents/SNPChip/STAG_Metadata_Template_v3.xlsm) for their samples to be analyzed.  A row is entered into the spreadsheet for each sample, and when finished, the spreadsheet is exported from Excel as a tab-delimeted file.  Naming the file is critical - the word *"metadata"* must be contained within the file name on disk.
 - The user uploads their sample metadata file (named something like *"affy_metadata.tabular"*) along with the necessary raw Affymetrix data files for their current run into the Galaxy CoralSNP environment using the *"Upload file"* tool within the *"Get Data"* section of the Galaxy tool panel.
 - The user selects the *"Queue genotype workflow"* tool from the *"Genotype Workflow"* section of the Galaxy tool panel, selectes the appropriate files as inputs, and executes the tool.  The user is now finished - the tool executes the enitre analysis for their samples.
