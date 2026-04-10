######################
#This is the R script to launch TransTag Shiny app to process raw sequencing files, and output the top enriched genomic sequences flanking transgene insertion sites.
#You can search/blast the most enriched flanking genomic sequence(s) in the genome to find the possible location of insertion site(s).
#author: Xiaolu Wei (xiaolu.wei@unt.edu)
######################

# Load libraries ----
library(shiny)
library(tidyverse)
library(dplyr)

#set upload input file size limit, default is set to 20 GB
options(shiny.maxRequestSize = 20 * 1000 * 1024^2)

remove_bold <-"#expr-container label {font-weight: 400;}"

# Define User Interface ----
ui <- fluidPage(
  titlePanel("TransTag"),
  
  sidebarLayout(
    ## Sidebar ----
    sidebarPanel(width = 4,
                 tags$style(remove_bold), 
                 tags$div(id = "expr-container", 
                   h4("Upload raw sequencing file"),
                   fileInput(inputId = "InputFile_raw",
                             label = "Accepted file formats are '.fastq', '.fq', '.fastq.gz' and '.fq.gz'",
                             accept = c(".gz", ".fastq", ".fq"),
                             multiple = FALSE),
                   h6(br()),
                   h4("Parameter"),
                   sliderInput(inputId = "cutoff", 
                               label = "Read length cutoff quantile", 
                               min = 0.2, max = 0.9, step= 0.05, value = 0.75),
                   textOutput(outputId = "size"),
                   h4(br()),
                   htmlOutput(outputId = "percentage"),
                 ) ## Closes custom style
      ), ## Closes Sidebar Panel
    ## Main Panel ----
    mainPanel(width = 8,
              tabsetPanel(
                
                ### Tab 1: Data Tables ----
                tabPanel("Output Table",
                         #h3("Output Table"),
                         column(width = 12,
                                h3("Top15 flanking genomic regions of insertion sites"),
                                "With the top flanking genomic region sequences, please use",
                                tags$a(href="https://genome.ucsc.edu/cgi-bin/hgBlat?hgsid=2483396421_8L6l8Cv6o6gPfhcNS6Ajh4JZT4JE&command=start", 
                                       "UCSC Genome Browser blast tool"), "to further locate the insertion site in the genome.",
                                tableOutput("out_simple_table"))
                ), ## Closes Tab 1
                
                ### Tab 2: Plots ----
                tabPanel("Summary Plot",
                         h3("Trimmed/processed chimeric read size distribution"),
                         plotOutput("out_plot")
                ) ## Closes Tab 2) ## Closes Tabset Panel
                
              ) ## Closes Tabset Panel
              
              ## Closing Brackets ----
    ) ## Closes Main Panel
  ) ## Closes Sidebar Layout
  
) ## Closes UI


# Define Server Logic ----
server <- function(input, output) {
  
  # Function Initiation ----
  k_count <- function(reads, k) {
    #only keep reads with length greater than kmer size
    filter_reads <- subset(reads, nchar(reads$V1) >= k)
    #extract kmers from remaining reads
    kmers <- data.frame(kmers = strtrim(filter_reads[,], k))
    #total number of kmers
    totalKmer <- nrow(kmers)
    
    #backbone sequence
    bg_seq <- strtrim("TCAAGAACTCCTGGACAAACCTCTGACCTGTGTGGAACAGAGTGGATATGGGTGTCTGAACAGATATTCACGTCTTTTGCAGATCAGAGGGCATTTCTGGTG", k)
    #kmer sequences that match the backbone sequence (with up to 5 mismatches)
    removed <- data.frame(kmers = unique(kmers[agrep(bg_seq, kmers[,], fixed = TRUE, max.distance = (all = 5)), ]))
    
    #remove kmers that match backbone sequences
    if (nrow(removed) > 0) {
      for (i in 1:nrow(removed))
      { kmers <- kmers %>% filter(.data[["kmers"]] != removed[i,])
      }
    }
    
    #number of kmers that do not match backbone sequence
    noBgKmer <- nrow(kmers)
    #percentage of backbone sequences in all kmers
    percentage <- 100*(1-noBgKmer/totalKmer)
    #number of kmers with backbone sequences
    number <- totalKmer-noBgKmer
      
    out <- kmers %>% group_by(kmers) %>% summarize(count = n()) %>% arrange(desc(count)) %>% slice(1:15)
    return(list(out, percentage, number))
  }
  
  ## Reactive expression to load data
  #input raw sequencing data file (fastq format)
  data <- reactive({
    req(input$InputFile_raw)
    #use system call to process the input file as following:
    #Extract chimeric reads that contain Tol2 sequences
    #Trim off Tn5 adapter sequences from 3' end
    #Trim off Tol2 sequences from 5' end
    inData <- as.data.frame(system(paste("zgrep TTTCACTTGAGTAAAATTTTTGAGTACTTTTTACACCTCTG", input$InputFile_raw$datapath, "| sed 's/CTGTCTCT.*$//g' | sed 's/^.*CTTTTTACACCTCTG//g'"), intern = TRUE))
    data.frame(V1=inData[,1])
  })
  
  ## Reactive expression to calculate read size cutoff
  size <- reactive({
    k <- as.numeric(input$cutoff)
    k_length <- quantile(nchar(data()$V1), k)
    paste0("Read size cutoff: ", as.character(k_length))
  })
  
  ## Reactive expression to calculate count
  count <- reactive({
    k <- as.numeric(input$cutoff)
    k_length <- quantile(nchar(data()$V1), k)
    k_count(data(), k_length)[[1]]
  })

  ## Reactive expression to calculate percentage/number of backbone kmers
  percentage <- reactive({
    k <- as.numeric(input$cutoff)
    k_length <- quantile(nchar(data()$V1), k)
    pt <- k_count(data(), k_length)[[2]]
    nn <- k_count(data(), k_length)[[3]]
    ifelse(pt > 30 & nn > 100, 
          paste0(paste0("<B>Warning</B>: Percentage of kmers with Tol2 vector sequences is ", sprintf("%0.1f%%", pt)), ", there may be Tol2-independent integration event(s)"),
          '')
    })
  
  ## Render text/table/plot output
  ## Text Outputs
  output$size <- renderText(size())
  output$percentage <- renderText(percentage())

  ## Table Outputs
  output$out_simple_table <- renderTable(count())

  ## Plot Output
  output$out_plot <- renderPlot({
    x <- nchar(data()$V1)
    hist(x, xlab = "Read size", ylab = "Counts", main = "")
  })
  
} ## Closes Server

# Run the Application ----
shinyApp(ui = ui, server = server)
