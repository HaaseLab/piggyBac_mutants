# *"Transposase N-terminal phosphorylation and asymmetric transposon ends inhibit piggyBac transposition in mammalian cells"*

### *Wentian Luo, Alison B. Hickman, Pavol Genzor, Rodolfo Ghirlando, Christopher M Furman, Anna Menshikh, Astrid D Haase, Fred Dyda, Matthew H Wilson*  





_Computational methods by: Pavol Genzor_

Accepted by NAR: November 2022
First deposited on Githhub: 06.24.2022

####
This document contains sample of computational methods associated with the manuscript. Please see the folders in this repository for sample data, mardowns/PDF of methods, and functions specific for this manuscript. Please contact Dr. Astrid Haase (astrid.haase@nih.gov) with any questions.   

Related code, sample data and functions are also available at github <https://github.com/HaaseLab/piggyBac_mutants> page.   


####
#### Table of Content

* Raw data processing
  * Read trimming and genomic alignments parameters
  
* ANALYSIS
  * Integration and target site duplication analysis
  * Loading paired-end (PE) .bam into R 
  * Identifying integration regions - "peaks"

* PLOTS
  * Chromosome peak distribution 
  * Genome annotation
  * Tiled integration map
  * Genome coverage
  
* Functions
  * filterBamPE()
  * findPeaksInPERegions()
  
