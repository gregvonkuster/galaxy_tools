Galaxy Tools maintained by [Greg Von Kuster](https://github.com/gregvonkuster)
==========================================

This repository contains tools that can be installed from the Galaxy Tool Shed [Galaxy Tool Shed](https://toolshed.g2.bx.psu.edu/) for use within [Galaxy](https://galaxyproject.org).

Highlights
----------
 * ChIP-seq - [Mahony lab](http://mahonylab.org)  at Penn State University

   * MultiGPS - a framework for analyzing collections of multi-condition ChIP-seq datasets and characterizing differential binding events between conditions.  MultiGPS encourages consistency in the reported binding event locations across conditions and provides accurate estimation of ChIP enrichment levels at each event.

 * Entomology - [Fleischer lab](https://ento.psu.edu/directory/sjf4)  at Penn State University

   * Temperature data for insect phenology model - data source tool for retrieving temperature data from the remote data source [Pestwatch](http://www.pestwatch.psu.edu/)
   * Insect phenology model - an agent-based stochastic model expressing stage-specific phenology and population dynamics for an insect species across geographic regions.
   * Extract date interval from insect phenology model data - extracts a date interval from the data produced by the Insect Phenology Model tool, providing a "zoomed in" view of the plots.

 * Epigenetics - [IDEAS pipeline](http://personal.psu.edu/yzz2/IDEAS/) from the [Zhang lab](https://yulili2000.wixsite.com/website)  at Penn State University

   * IDEAS Preprocessor - maps a list of epigenetic datasets to a common genomic coordinate in a selected assembly, producing datasets for use as input to IDEAS.
   * IDEAS - an Integrateive and Discriminitive Epigenome Annotation System that identifies de novo regulatory functions from epigenetic data in multiple cell types jointly.
   * IDEAS Genome Tracks - creates [UCSC Genome Browser Track Hubs](https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html) for vizualizing IDEAS outputs.

 * PlantTribes - [PlantTribes pipelines](https://github.com/dePamphilis/PlantTribes/tree/master/pipelines) from the [DePamphillis lab](http://cwd.huck.psu.edu/) at Penn State University.

   * Load PlantTribes Scaffold - analyzes scaffolds installed into Galaxy by the PlantTribes Scaffolds Downloader data manager tool and inserts information about them into the Galaxy PlantTribes database for querying and additional analysis.
   * Update PlantTribes Scaffold - adds a new genome to a scaffold installed into Galaxy by the PlantTribes Scaffolds Downloader data manager tool.
   * AssemblyPostProcessor - post-processes de novo assembled transcripts into putative coding sequences and their corresponding amino acid translations and optionally assigns transcripts to circumscribed gene families (orthogroups).
   * GeneFamilyClassifier - classifies gene coding sequences either produced by the AssemblyPostProcessor tool or from an external source into pre-computed orthologous gene family clusters (orthogroups) of a PlantTribes scaffold.
   * GeneFamilyIntegrator - integrates PlantTribes scaffold orthogroup backbone gene models with gene coding sequences classified into the scaffold by the GeneFamilyClassifier tool.
   * GeneFamilyAligner -  estimates protein and codon multiple sequence alignments of integrated orthologous gene family fasta files produced by the GeneFamilyIntegrator tool.
   * GeneFamilyPhylogenyBuilder - performs gene family phylogenetic inference of multiple sequence alignments produced by the GeneFamilyAligner tool.
   * KaKsAnalysis - estimates paralogous and orthologous pairwise synonymous (Ks) and non-synonymous (Ka) substitution rates for a set of gene coding sequences either produced by the AssemblyPostProcessor tool or from an external source.
   * KsDistribution - uses the analysis results produced by the KaKsAnalysis tool to plot the distribution of synonymous substitution (Ks) rates and fit the estimated significant normal mixtures component(s) onto the distribution.
