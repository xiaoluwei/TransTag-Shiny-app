# TransTag

This repository contains scripts for the alignment-free Shiny app analysis described in: 

**TransTag enables simple and efficient transgene mapping in zebrafish via tagmentation.** <br/>
Fanju W. Meng, Paige Schneider, Xiaolu Wei, Krishan Ariyasiri, Marnie E. Halpern, Patrick J. Murphy. 2025. [Cell Reports Methods](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(25)00126-2)


## Script 

**TransTag_alignmentFree.ShinyApp.R**

This is the R script to launch Shiny app interface to process the raw sequencing reads fastq.gz file (first read R1 file for paired-end reads) from the TransTag library. Briefly, the script extracts chimeric reads that contain Tol2 sequences, trim offs Tn5 adapter and Tol2 sequences to extract the genomic sequences flanking transgene insertion sites, and displays the top fifteen most abundant flanking region sequences. The enriched flanking sequences can then be mapped to the zebrafish genome using any standard online tool, such as UCSC BLAT (https://genome.ucsc.edu) or NCBI BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi), allowing users to assign genomic coordinates of transgenes. <br/>

- Example input file: ```Example_subsample_to_1millionReads_R1.fastq.gz``` <br/>

- Users can use the online version of TransTag Shiny app interface (for smaller size input files (<1 GB), the shinyapps.io can only process input files smaller than 1GB). <br/>

  - a. Open alignment-free R Shiny App website in a web browser: https://menglab.shinyapps.io/transtag_alignmentfree/. <br/>
  - b. Upload raw sequencing data from the TransTag library for processing (first read R1 file for paired-end reads). It will show the top fifteen abundant genomic sequences flanking Tol2 insertion site in the “Output Table” tab, and size distribution of the chimeric reads with flanking genomic sequences in the “Summary Plot” tab. You may choose to change the read length cutoff quantile between 0.2 to 0.9 based on the size distribution. The default value for read length cutoff quantile is set at 0.75. <br/>
  - c. By changing the read length cutoff quantile between 0.2 to 0.9, the TransTag Shiny app will display flanking genomic regions with different lengths and count numbers accordingly. <br/>
  - d. To identify the genomic coordinates of transgene insertion site(s), use any standard online tool, such as UCSC BLAT (https://genome.ucsc.edu) or NCBI BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi). <br/>
  
- Alternatively, users can run TransTag_alignmentFree.ShinyApp.R on their own computers to launch the Shiny app interface (for any size input files). <br/>

  - a. Download and install R and RStudio on local computer following instruction on https://posit.co/download/rstudio-desktop/. <br/>
  - b. Download alignment-free R Shiny app script “TransTag_alignmentFree.ShinyApp.R”. <br/>
  - c. Double click the downloaded R shiny script to open it in RStudio. <br/>
  - d. Install required packages. Type the following code in RStudio Console and enter to run: <br/>
  	 ```
	 install.packages("shiny")
	 install.packages("tidyverse")
	 install.packages("dplyr")
	 ```
  - e. Click “Run App” to launch the Shiny app. The Shiny app interface will show up in a new window. Please follow the same steps as the using the online version of TransTag Shiny app interface section to conduct analysis. <br/>


## Notes

1. Based on the assembled Tn5 used in the library preparation step, R1 reads file for the pair-end sequencing reads will have the Tol2 repeat sequence. For alignment-free analysis, R1 reads file would be the input file for the Shiny app interface. <br/>
2. If "Maximum upload size exceeded" error occurs, modify `options(shiny.maxRequestSize = )` in the TransTag_alignmentFree.ShinyApp.R script to increase the input file size limit. The default size limit is 20 GB. <br/>
3. Please refer to the "TransTag_Shinyapp_analysis_tutorial" in this repository for more detailed step-by-step tutorial for data analysis.
