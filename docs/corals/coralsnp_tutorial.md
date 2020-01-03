# CoralSNP - a Galaxy analysis environment that includes a 30k SNP genotyping array for Acropora corals and their dinoflagellate symbionts

## Introduction

This document provides information for using the [CoralSNP Galaxy environment](https://coralsnp.science.psu.edu/galaxy)
which is based on the [Galaxy workbench](https://galaxyproject.org/).  A general understanding of Galaxy is required, so please spend some time with the [introductory tutorials](https://training.galaxyproject.org/training-material/topics/introduction) if you are not yet familiar with Galaxy.

The Galaxy CoralSNP environment enables streamlined analysis of the coral SNPchip available from Fisher Scientific to ultimately provide the user with a genet id, converted raw genotyped data, sample relatedness and hybrid status.

The process is straightforward.  Each of the following steps will be discussed in detail in the following sections.

 - A sample metadata file is created by the user from an Excel spreadsheet template for their samples to be analyzed.  A row is entered into the spreadsheet for each sample, and when finished, the spreadsheet is exported from Excel and saved to disk as a tab-delimeted file.  Naming the file is critical - the word *"metadata"* must be contained within the file name on disk.
 - The user logs into the [CoralSNP Galaxy](https://coralsnp.science.psu.edu/galaxy) environment and creates a new, empty history, ideally naming it in a way that associates it with the run being analyzed.
 - The user uploads their sample metadata file (named something like *"affy_metadata.tabular"*) along with the necessary raw Affymetrix data files for the run being analyzed into the Galaxy CoralSNP environment using the *"Upload file"* tool within the *"Get Data"* section of the Galaxy tool panel.
 - The user selects the *"Queue genotype workflow"* tool from the *"Genotype Workflow"* section of the Galaxy tool panel, selects the appropriate files as inputs, and executes the tool.  The tool executes the entire analysis for the samples and the user can view the results in the Galaxy history when the analysis is finished.

## Creating the Sample Metadata File

The metadata file for the run describes the samples being analyzed by providing important information about them.  The Baums' Lab website provides an [Excel spreadsheet template](http://baumslab.org/documents/SNPChip/STAG_Metadata_Template_v3.xlsm) that can be downloaded and used for each sample run.  Some of the data is optional.

 - **user_specimen_id**		(required)
 - **field_call**		(optional)
 - **bcoral_genet_id**		(optional) - the Baums' lab coral genet id, deprecated but remains for backward compability with earlier analyses
 - **bsym_genet_id**		(optional) - the Baums' lab symbiont id, deprecated but remains for backward compability with earlier analyses
 - **reef**			(required) - the name of the reef from which the samples were collected
 - **region**			(required) - the geographic region in which the reef is located
 - **latitude**			(required) - the latitude where the sample was collected
 - **longitude**		(required) - the longitude where the sample was colected
 - **geographic_origin**	(optional) - 
 - **colony_location**
 - **depth**
 - **disease_resist**
 - **bleach_resist**
 - **mortality**
 - **tle**
 - **spawning**
 - **collector_last_name**
 - **collector_first_name**
 - **org**
 - **collection_date**
 - **contact_email**
 - **seq_facility**
 - **array_version**
 - **public**
 - **public_after_date**
 - **sperm_motility**
 - **healing_time**
 - **dna_extraction_method**
 - **dna_concentration**
 - **registry_id**
 - **result_folder_name**
 - plate_barcode**


It is crucial to ensure that the information in this sample metadata file is correct.  
