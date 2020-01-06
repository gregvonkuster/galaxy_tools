# CoralSNP - a Galaxy analysis environment that includes a 30k SNP genotyping array for Acropora corals and their dinoflagellate symbionts

## Introduction

This document provides information for using the [Galaxy CoralSNP environment](https://coralsnp.science.psu.edu/galaxy)
which is based on the [Galaxy workbench](https://galaxyproject.org/).  A general understanding of Galaxy is required, so please spend some time with the [introductory tutorials](https://training.galaxyproject.org/training-material/topics/introduction) if you are not yet familiar with Galaxy.

The Galaxy CoralSNP environment enables streamlined analysis of the coral SNPchip available from Fisher Scientific to ultimately provide the user with a genet id, converted raw genotyped data, sample relatedness and hybrid status.

The process is straightforward.  Each of the following steps will be discussed in detail in the following sections of this document.

 - A sample metadata file is created by the user from an Excel spreadsheet template for their samples to be analyzed.  A row is entered into the spreadsheet for each sample, and when finished, the spreadsheet is exported from Excel and saved to disk as a tab-delimeted file.
 - The user logs into the Galaxy CoralSNP environment and creates a new, empty history, ideally naming it in a way that associates it with the run being analyzed.
 - The user uploads their sample metadata file along with the necessary raw Affymetrix data files into the Galaxy CoralSNP environment using the *"Upload file"* tool within the *"Get Data"* section of the Galaxy CoralSNP tool panel.
 - The user selects the *"Queue genotype workflow"* tool from the *"Genotype Workflow"* section of the Galaxy CoralSNP tool panel, selects the appropriate files as inputs, and executes the tool.  The tool executes the entire analysis for the samples and the user can view the results in the Galaxy history when the analysis is finished.

## Create the Sample Metadata File

The metadata file for the run describes the samples being analyzed by providing important information about them.  The Baums' Lab website provides an [Excel spreadsheet template](http://baumslab.org/documents/SNPChip/STAG_Metadata_Template_v3.xlsm) that can be downloaded and used for each sample run.  Some of the data is optional - here are some details about the columns in the spreadsheet template.

 - **user_specimen_id** (required) - user-specific identifier for each sample
 - **field_call** (optional)
 - **bcoral_genet_id** (optional) - the Baums' lab coral genet id, deprecated but remains for backward compability with earlier analyses
 - **bsym_genet_id** (optional) - the Baums' lab symbiont genet id, deprecated but remains for backward compability with earlier analyses
 - **reef** (required) - the name of the reef from which the samples were collected
 - **region** (required) - the geographic region in which the reef is located
 - **latitude** (required) - the latitude (in decimal degrees) where the sample was collected
 - **longitude** (required) - the longitude (in decimal degrees) where the sample was colected
 - **geographic_origin** (optional) - where the sample was taken (must be either *"colony"* or *"reef"*), defaults to *"reef*"
 - **colony_location** (optional)
 - **depth** (optional) - depth (in meters) from surface where the samples was taken, must be decimal value
 - **disease_resist** (optional)
 - **bleach_resist** (optional)
 - **mortality** (optional)
 - **tle** (optional)
 - **spawning** (optional)
 - **collector_last_name** (required) - the last name of the collector
 - **collector_first_name** (required) - the first name of the collector
 - **org** (required) - the organization for which the collector is working
 - **collection_date** (required) - the date (format yyyy-mm-dd) when the sample was collected
 - **contact_email** (required) - the collector's email address
 - **seq_facility** (required) - the facility sequencing the samples
 - **array_version** (optional)
 - **public** (optional) - whether the information about the samples is public, defaults to *"Yes"*
 - **public_after_date** (optional) - the date at which information about the samples can be made public, defaults to the current day
 - **sperm_motility** (required)
 - **healing_time** (required)
 - **dna_extraction_method** (optional) - method used to extract the sample's dna
 - **dna_concentration** (optional) - dna concentration in ug
 - **registry_id** (optional)
 - **result_folder_name** (optional)
 - **plate_barcode** (optional)


It is crucial that the information in this sample metadata file is correct since the analysis pipeline will store portions of it in both current and future analyses for comparison and other uses.  The Excel spreadsheet doesn't validate the information, so care must be taken when entering the data.

Information about all previous runs is available in the [CoralSNP Reports](https://coralsnp.science.psu.edu/reports).  Users should browse this information before creating their sample metadata file to help ensure that the information is correct.  For example, the collector's name, organization and email, reef name, region, etc should exactly match those entered for previous samples if taken by the same collector, from the same reef, etc.

When the Excel spreadsheet is complete. export the information (tab-separated format) into a disk file.  Naming the file is critical - the word *"metadata"* must be contained within the file name on disk.  For example, if analyzing Plate 2 data, the file could be named *"plate2_metadata.tabular"*, but any file name that includes the string *"metadata*" is fine (e.g., *"affy_metadata.tabular"*).

## Upload the Sample Data to Galaxy CoralSNP for Analysis

You'll have to create an account upon your initial visit to the [Galaxy CoralSNP](https://coralsnp.science.psu.edu/galaxy) environment.  General information about registering and logging into Galaxy is available [here](https://galaxyproject.org/support/account), and can be leveraged for creating a new account in the Galaxy CoralSNP environment for users that have not yet done so.  If logging in for the first time, a new, empty Galaxy history will be created for you.  If you are logging into an existing account, you should create a new, empty history for the samples you're planning to analyze.  Name the history in a way that associates it with your current samples.

It is imperative that a new analysis is performed within a Galaxy history that contains only the dataset items for the run.  Analyses should not be performed within a history that contains any items except for the input datasets for the current run.  The analysis pipeline inspects the current Galaxy history and specifies named items as designated inputs to certain tools within the pipeline, so multiple items with the same name will result in pipeline errors.  In addition, all items within the history must be in the *"ok"* (i.e., *"green"*) state.  History items in a *"queued"*, *"running"* or *"error"* state will result in analysis errors and cause problems for others attempting to perform analyses.  Problematic items can be deleted from a history by clicking the *X* icon in the item, allowing for histories to be "cleaned up" in preparation for performing an analyis.  Wenever possible, items being deleted from a history should be *"permamently deleted"* as this will free up disk space on the Galaxy CoralSNP server.

You can now upload your data for analysis.  A general tutorial for uplaoding data to Galaxy is avaialble [here](https://training.galaxyproject.org/training-material/topics/galaxy-data-manipulation/tutorials/get-data/slides.html#1).

Here is a view of the Galaxy upload form that shows all of the files for what was called *"Plate 2"* of the samples.  Notice that the data type (e.g., *"tabular"*. *"csv"*, *"txt"*, etc) has been selected for each file.  Although not required (Galaxy will auto-detect file formats), this is ideal since it will descrease the time needed to upload all of the data files.  Also notice the names of the files.  With the exception of the *"affy_metadata.tabular"* file which is named by the user, all of the file names are produced by the Affymetrix sequencing process for the samples.  All of these file names are important - as discussed previously, the analysis pipeline specifies these named files as designated inputs to certain tools within the pipeline.

![Raw Affymetrix sample data and samples metadata file for Plate 2](upload_data.png)

When all of the files have been chosen and the data formats specified, clicking the *"Start"* button will initiate the upload.

## The *"Queue genotype workflow"* Tool

After uploading the files, select the *"Queue genotype workflow"* tool from the *"Genotype Workflow"* section of the Galaxy CoralSNP tool panel.  Here is a view of the tool form and the Galaxy history containing the uploaded Plate 2 sample files.  Notice that the Galaxy history contains only the uploaded files, and that all items are in the *"ok"* (i.e., *"green"*) state as discussed previously.  This is essential in order to ensure a successful analysis.

![Specify appropritate inputs](queue_genotype_workflow.png)

We mentioned previously that the names of the files are important, and the tool form demonstrates this.  Notice the help text below each input selection (e.g., *"The word 'metadata' must be in the file name."*).

Executing this tool invokes the entire analysis pipeline which consists of a complex set processes and components.
