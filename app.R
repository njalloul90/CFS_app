#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)

library(DT)
library(data.table)
library(ggbio)
library(plyr)
library(dplyr)

library(stringr)
library(karyoploteR)

library(gtools)
# library(BiocManager)
# options(repos = BiocManager::repositories())
library(BiocManager)
library(BSgenome.Hsapiens.UCSC.hg38) 
library(BSgenome.Hsapiens.UCSC.hg19) 
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(regioneR)
library(Rsamtools)
library(patchwork)

library("ggplot2") # for the plot
library("ggrepel") # for spreading text labels on the plot
library("scales") # for axis labels notation

library(qualV)

bases <- c("C", "G", "A", "T")

# Define UI for application that draws a histogram
ui <- fluidPage(
    shinyjs::useShinyjs(),
    # Application title
    h3("Common Fragile Sites in the Human Genome"),
    h6("The interactive display features Common Fragile Sites (CFSs) in the reference human genome."),
    h6("It allows you to explore the DNA sequences that forms these regions and compare them to other non-CFSs in the same chromosome."),
    HTML('<hr style="color: black;">'),
    
    # Panel 1 : Ref Seq, Chr, MotifL
    fluidRow(column(3,
                    # reference genome selection
                    selectInput("refGenome", "Reference Genome:",
                                c("GRCh37 hg19" = "hg19",
                                  "GRCh38 hg38" = "hg38",
                                  "T2T-CHM13v2.0" = "T2T"))
                    ),
             column(3,
                    # chromosome selection
                    selectInput("chr", "Chromosome:",
                                c("Chr1" = "chr1", "Chr2" = "chr2", "Chr3" = "chr3", "Chr4" = "chr4",
                                  "Chr5" = "chr5", "Chr6" = "chr6", "Chr7" = "chr7", "Chr8" = "chr8",
                                  "Chr9" = "chr9", "Chr10" = "chr10", "Chr11" = "chr11", "Chr12" = "chr12",
                                  "Chr13" = "chr13", "Chr14" = "chr14", "Chr15" = "chr15", "Chr16" = "chr16",
                                  "Chr17" = "chr17", "Chr18" = "chr18", "Chr19" = "chr19", "Chr20" = "chr20",
                                  "Chr21" = "chr21", "Chr22" = "chr22", "ChrX" = "chrX", "ChrY" = "chrY")
                                )
                    ),
             column(3,
                    # Motif Length selection
                    selectInput("motifL", "Motif Length:",
                                c("L = 3" = "3", "L = 4" = "4", "L = 5" = "5",
                                  "L = 6" = "6", "L = 7" = "7", "L = 8" = "8"))
                    ),
             column(3,
                    # text output n possible permutations of motif length L
                    textOutput("nMotifs")
                    )
             ),
    
    HTML('<hr style="color: black;">'),
    
    # Panel 2 : Chr Coordinates , Fragile Sites for selected chr
    fluidRow(
        column(3,h5("Selected Chromosome"),tableOutput('tableChrCoord')),
        column(6,DT::dataTableOutput('tableFragileSites'))
                    ),
    
    HTML('<hr style="color: black;">'),
    
    # Panel 3 : Either upload table or choose selected input for AT_GC content calculations
    fluidRow(column(3,
                    # Upload table for AT_GC content
                    fileInput("userInputTable", "Choose CSV/TXT File",
                              multiple = TRUE,
                              accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
                    h6("N.B. select appropriate chr and reference sequence"),
                    h6("OR"),
                    checkboxInput("chooseSelectedTable", "Choose Selected Input", FALSE),
                    h6("N.B. unselect to enable user file upload.")
                    ),
             column(6,
                    # Display AT-GC Content Table
                    h5("GC-AT Content"),
                    h6("Requireeed Headers: chr, StartPos, EndPos, Label"),
                    DT::dataTableOutput("AT_Table")
                    )),
    
    HTML('<hr style="color: black;">'),
    # Panel 4: plot nucleotide frequency across regions
    fluidRow(column(3,
                    h5("Nucleotide Frequency"),
                    numericInput("nuclFreqInterval", "Nucleotide frequency per bp:", 10000, min = 1000, max = 1000000,step = 1000),
                    h6("N.B. minimum interval 1Kb, max interval 1Mb"),
                    actionButton("button_karyo", "Plot karyotype")
                    ),
             column(6,
                    plotOutput("plot_karyo",height = 300)
                    )
             ),
    HTML('<hr style="color: black;">'),
    # Panel 5 Motif counts
    fluidRow(column(3,
                    h5("Motifs"),
                    h6("Motifs are generated from {'A','T','C','G'} combinations."),
                    h6("Get Motif Counts: calculates the number of times a motif exists in the defined reegions."),
                    h6("Get Motif Positions: finds the exact chromosomal coordinates of the motifs."),
                    HTML('<hr style="color: black;">'),
                    radioButtons("motifOverlap", "Overlap Motifs", 
                                 list("Yes","No"), inline = TRUE, selected = "Yes"),
                    h6("Overlap Motifs:"),
                    h6("e.g. in seq 'ATTATATTAT'; motif 'TAT' count = 3 if motif overlap is chosen."),
                    h6("if no overlap is chosen; motif 'TAT' count = 2."),
                    actionButton("button_MotifCounts", "Get Motif Counts"),
                    HTML('<hr style="color: white;">'),
                    #actionButton("button_MotifPositions", "Get Motif Positions"),
                    HTML('<hr style="color: white;">'),
                    uiOutput("regionToAdd")
                    ),
             column(6,
                    DT::dataTableOutput("MotifCounts_Table",width = 900)
                    #DT::dataTableOutput("MotifPositions_Table",width = 900)
                    )),
    HTML('<hr style="color: black;">'),
    fluidRow(column(3,
                    h5("Upload Motif Counts Table"),
                    fileInput("motifCountsTable", "Choose CSV File",
                              multiple = TRUE,
                              accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
                    h5("Select Inputs From Column Headers"),
                    selectInput('MotifNames','Motifs',""),
                    selectInput('Counts_ControlRegion', 'Count in Control Region', ""),
                    selectInput('norm_Counts_ControlRegion', 'Norm Count in Control Region', ""),
                    selectInput('Count_Repeats_ControlRegion', 'Count Repeats in Control Region', ""),
                    selectInput('norm_Count_Repeats_ControlRegion', 'Norm Count Repeats in Control Region', ""),
                    selectInput('countROI', 'Counts in ROI', ""),
                    selectInput('norm_countROI', 'Norm Counts in ROI', ""),
                    selectInput('repeat_countROI', 'Count Repeats in ROI', ""),
                    selectInput('norm_repeat_countROI', 'Norm Count Repeats in ROI', ""),
                    numericInput("sizeControl", "Size of Control Region:", value= 1000, min = 1, max = 100000000),
                    numericInput("sizeROI", "Size of ROI:", value = 1000, min = 1, max = 100000000),
                    HTML('<hr style="color: white;">'),
                    actionButton("Stats_MotifsPresentHigherRate", "Get Motifs Present at Higher Rate in ROI"),
                    HTML('<hr style="color: white;">'),
                    actionButton("Stats_MotifsRepeatedHigherRate", "Get Motifs Repeated at Higher Rate in ROI")),
             column(9,DT::dataTableOutput("MotifCountsTableUpload"))),
    HTML('<hr style="color: black;">'),
    fluidRow(column(3,
                    h5("Motifs Present at a Higher Rate in ROI")),
             column(6,
                    DT::dataTableOutput("Motifs_Present_Higher_Rate"))
             ),
    HTML('<hr style="color: black;">'),
    fluidRow(column(3,
                    h5("Motifs Repeated at a Higher Rate in ROI")),
             column(6,
                    DT::dataTableOutput("Motifs_Repeated_Higher_Rate"))
    ),
    HTML('<hr style="color: black;">')
    
)

server <- function(input, output,session) {
    
    # nMotifs: text output for n possible permutations of motif length L
    output$nMotifs <- renderText({
        paste("Possible permutations of motif length L:", 4^strtoi(input$motifL))
    })
    
    # tableChrCoord: get selected reference genome and display chr coordinates accordingly
    output$tableChrCoord <- renderTable({
        # assign reference genome
        if (input$refGenome == "hg19") {
            refSeq <- "hg19"
            BSGENOME_SELECTED = BSgenome.Hsapiens.UCSC.hg19
            # assign chromosome from selection
            chrLengths = seqlengths(BSGENOME_SELECTED)
            Chr = input$chr
            Start = 1
            End = chrLengths[[input$chr]]
        }
        if (input$refGenome == "hg38") {
            refSeq <- "hg38"
            BSGENOME_SELECTED = BSgenome.Hsapiens.UCSC.hg38
            # assign chromosome from selection
            chrLengths = seqlengths(BSGENOME_SELECTED)
            Chr = input$chr
            Start = 1
            End = chrLengths[[input$chr]]
        }
      if (input$refGenome == "T2T") {
        refSeq <- "T2T"
        fasta_file = paste0(getwd(),"/Data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna")
        fa = FaFile(fasta_file) 
        BSGENOME_SELECTED = as(seqinfo(fa), "GRanges")
        chrLengths = BSGENOME_SELECTED@seqinfo@seqlengths
        Chr = input$chr
        Start = 1
        End = chrLengths[strtoi(str_remove(Chr, "chr"), base=0L)]
        
        BSGENOME_SELECTED = BSGENOME_SELECTED
        
      }
        
        tableChrCoord = data.frame(refSeq,Chr,Start,End)
    })
    
    # tableFragileSites: get Fragile sites for the selected chromosome and reference
    output$tableFragileSites <- DT::renderDataTable({
        # load FS_COORD Table
        FS_COORD = "FS_COORD.txt"
        FS_COORD <- read.table(FS_COORD, header=TRUE, sep = "\t")
        # filter FS based on selected refseq and chr
        # FS_COORD_REFSEQ_CHR = FS_COORD[FS_COORD$Chromosome == input$chr,]
        # FS_COORD_REFSEQ_CHR = FS_COORD[FS_COORD$REF == input$refGenome,]
        FS_COORD_REFSEQ_CHR = subset(FS_COORD, Chromosome == input$chr & REF == input$refGenome)
    })
    
    # AT_Table: output GC-AT content based on userInputTable or on chooseSelectedTable
    output$AT_Table <- DT::renderDataTable(extensions = 'Buttons', options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
    ),{
        # if checkbox is selected, disable upload button
        observeEvent(input$chooseSelectedTable, {
            if (input$chooseSelectedTable) {
                shinyjs::disable("userInputTable")
            }
            if(!input$chooseSelectedTable) {
                shinyjs::enable("userInputTable")
            }
        })
        
        # if checkbox is selected, use coordinates from Panel 1
        if (input$chooseSelectedTable) {
            
            FS_COORD = "FS_COORD.txt"
            FS_COORD <- read.table(FS_COORD, header=TRUE, sep = "\t")
            FS_COORD_REFSEQ_CHR = subset(FS_COORD, Chromosome == input$chr & REF == input$refGenome)
            FS_COORD_REFSEQ_CHR <- FS_COORD_REFSEQ_CHR[order(FS_COORD_REFSEQ_CHR$StartPos),]
            
            # get seqLengths from selected refGenome
            if (input$refGenome == "hg19") {
                refSeq <- "hg19"
                BSGENOME_SELECTED = BSgenome.Hsapiens.UCSC.hg19
                # assign chromosome from selection
                chrLengths = seqlengths(BSGENOME_SELECTED)
                Chr = input$chr
                Start = 1
                End = chrLengths[[input$chr]]
                # get chromosome sequence
                chrseq = getSeq(BSGENOME_SELECTED, Chr,
                                start = Start, end = End)
            }
            if (input$refGenome == "hg38") {
                refSeq <- "hg38"
                BSGENOME_SELECTED = BSgenome.Hsapiens.UCSC.hg38
                # assign chromosome from selection
                chrLengths = seqlengths(BSGENOME_SELECTED)
                Chr = input$chr
                Start = 1
                End = chrLengths[[input$chr]]
                # get chromosome sequence
                chrseq = getSeq(BSGENOME_SELECTED, Chr,
                                start = Start, end = End)
            }
            if (input$refGenome == "T2T") {
              # T2T genome
              refSeq <- "T2T"
              fasta_file = paste0(getwd(),"/Data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna")
              fa = FaFile(fasta_file) 
              BSGENOME_SELECTED = as(seqinfo(fa), "GRanges")
              chrLengths = BSGENOME_SELECTED@seqinfo@seqlengths
              Chr = input$chr
              Start = 1
              End = chrLengths[strtoi(str_remove(Chr, "chr"), base=0L)]
              Chr_name_t2t = BSGENOME_SELECTED@seqinfo@seqnames
              Chr_name_t2t = Chr_name_t2t[strtoi(str_remove(Chr, "chr"), base=0L)]
              # get chromosome sequence
              chrseq = getSeq(fa, BSGENOME_SELECTED[strtoi(str_remove(Chr, "chr"), base=0L)],
                              start = Start, end = End)
            }
            
            
            # initialize DF
            df = data.frame(matrix(ncol = 4, nrow = 0))
            colnames(df) <- c("chr","StartPos","EndPos","Label")
            df[1,1] = Chr
            df[1,2] = 1
            df[1,3] = FS_COORD_REFSEQ_CHR$StartPos[1] - 1
            df[1,4] = "control"
            L = nrow(FS_COORD_REFSEQ_CHR)*2+1
            
            df[L,1] = Chr
            df[L,2] = FS_COORD_REFSEQ_CHR$EndPos[nrow(FS_COORD_REFSEQ_CHR)] +1
            df[L,3] = End
            df[L,4] = "control"
            
            # fill in remainder of df
            for (i in seq(2,L-1,by=2)) {
                df[i,1] = Chr
                df[i,2] = FS_COORD_REFSEQ_CHR$StartPos[i/2]
                df[i,3] = FS_COORD_REFSEQ_CHR$EndPos[i/2]
                df[i,4] = toString(FS_COORD_REFSEQ_CHR$FS[i/2])
            }
            for (i in seq(3,L-2,by=2)) {
                df[i,1] = Chr
                df[i,4] = "control"
                
                df[i,2] = df[i-1,3]+1
                df[i,3] = df[i+1,2]-1
            }
            df <- subset(df, EndPos>StartPos) 
            df$GC_Content = rep(0,nrow(df))
            df$AT_Content = rep(0,nrow(df))
            # calculate GC and AT content over each interval
            withProgress(message = 'Processing', value = 0, {
            for (row in 1:nrow(df)) {
                # Increment the progress bar, and update the detail text.
                incProgress(1/i, detail = paste("Calculating GC/AT content ", i))
                
                start_i = df[row,"StartPos"]
                end_i = df[row,"EndPos"]
                subseq = substring(chrseq,start_i,end_i)
                num_g <- str_count(subseq, "G")
                num_c <- str_count(subseq, "C")
                num_a <- str_count(subseq, "A")
                num_t <- str_count(subseq, "T")
                
                num_N <- str_count(subseq, "N")
                gc_content <- (num_g + num_c) / (str_length(subseq)-num_N) * 100
                at_content <- (num_a + num_t) / (str_length(subseq)-num_N) * 100
                df[row,"GC_Content"] = gc_content
                df[row,"AT_Content"] = at_content
            }
            })
            
            df
        }
        
        # if user inputs a file
        if (!input$chooseSelectedTable) {
            # get user input file
            userInputFile <- input$userInputTable
            ext <- tools::file_ext(userInputFile$datapath)
            req(userInputFile)

            # check file extenstion
            if (ext=="txt") {
                # read text file
                df <- read.table(userInputFile$datapath, header=TRUE, sep = "\t")
            }
            if (ext=="csv") {
                df <- read.table(userInputFile$datapath, header=TRUE, sep = ",")
            }
            if (ext!="txt" && ext!="csv") {

                df = data.frame(matrix(ncol = 6, nrow = 0))
            }
            
            # if not empty df, check colnames
            if (nrow(df) != 0) {
                requiredColumns <- c("chr","StartPos","EndPos","Label")
                inputColNames = colnames(df)
                # get seqLengths from selected refGenome
                if (input$refGenome == "hg19") {
                    refSeq <- "hg19"
                    BSGENOME_SELECTED = BSgenome.Hsapiens.UCSC.hg19
                    # assign chromosome from selection
                    chrLengths = seqlengths(BSGENOME_SELECTED)
                    Chr = input$chr
                    Start = 1
                    End = chrLengths[[input$chr]]
                    # get chromosome sequence
                    chrseq = getSeq(BSGENOME_SELECTED, Chr,
                                    start = Start, end = End)
                }
                if (input$refGenome == "hg38") {
                    refSeq <- "hg38"
                    BSGENOME_SELECTED = BSgenome.Hsapiens.UCSC.hg38
                    # assign chromosome from selection
                    chrLengths = seqlengths(BSGENOME_SELECTED)
                    Chr = input$chr
                    Start = 1
                    End = chrLengths[[input$chr]]
                    # get chromosome sequence
                    chrseq = getSeq(BSGENOME_SELECTED, Chr,
                                    start = Start, end = End)
                }
                
                if (input$refGenome == "T2T") {
                  # T2T genome
                  refSeq <- "T2T"
                  fasta_file = paste0(getwd(),"/Data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna")
                  fa = FaFile(fasta_file) 
                  BSGENOME_SELECTED = as(seqinfo(fa), "GRanges")
                  chrLengths = BSGENOME_SELECTED@seqinfo@seqlengths
                  Chr = input$chr
                  Start = 1
                  End = chrLengths[strtoi(str_remove(Chr, "chr"), base=0L)]
                  Chr_name_t2t = BSGENOME_SELECTED@seqinfo@seqnames
                  Chr_name_t2t = Chr_name_t2t[strtoi(str_remove(Chr, "chr"), base=0L)]
                  # get chromosome sequence
                  chrseq = getSeq(fa, BSGENOME_SELECTED[strtoi(str_remove(Chr, "chr"), base=0L)],
                                  start = Start, end = End)
                }
                
                if (length(intersect(requiredColumns, inputColNames))==length(requiredColumns)) {
                    # calculate GC-AT content
                    withProgress(message = 'Processing', value = 0, {
                        
                        for (row in 1:nrow(df)) {
                            # Increment the progress bar, and update the detail text.
                            incProgress(1/row, detail = paste("Calculating GC/AT content ", row))
                            
                            start_i = df[row,"StartPos"]
                            end_i = df[row,"EndPos"]
                            subseq = substring(chrseq,start_i,end_i)
                            num_g <- str_count(subseq, "G")
                            num_c <- str_count(subseq, "C")
                            num_a <- str_count(subseq, "A")
                            num_t <- str_count(subseq, "T")
                            
                            num_N <- str_count(subseq, "N")
                            gc_content <- (num_g + num_c) / (str_length(subseq)-num_N) * 100
                            at_content <- (num_a + num_t) / (str_length(subseq)-num_N) * 100
                            df[row,"GC_Content"] = gc_content
                            df[row,"AT_Content"] = at_content
                        }
                    })
                }
            }
            df
        }
        AT_Table = df
        write.table(AT_Table, 'AT_Table.txt', append = FALSE, sep = "\t", dec = ".",
                    row.names = FALSE, col.names = TRUE)
        #df_AT_Table <<- AT_Table
        AT_Table
    })
    
    
    # plot_karyo: plot karyotype based on nuclFreqInterval, and when button_karyo is pressed
    input_nuclFreqInterval <- eventReactive(input$button_karyo, {
        input$nuclFreqInterval
    })
    
    output$plot_karyo <- renderPlot({
      # make sure input$nuclFreqInterval is between 1000 and 1000000
      nuclFreq = input_nuclFreqInterval()
      if (!between(nuclFreq, 1000, 1000000)) {
        showNotification("Nucleotide frequency interval not between 1000 and 1000000 bp!")
      }
      
      if (between(nuclFreq, 1000, 1000000)) {
        # calculate frequency of nucleotides per interval
        if (input$refGenome=="hg19") {
          refSeq <- "hg19"
          BSGENOME_SELECTED = BSgenome.Hsapiens.UCSC.hg19
          # assign chromosome from selection
          chrLengths = seqlengths(BSGENOME_SELECTED)
          Chr = input$chr
          Start = 1
          End = chrLengths[[input$chr]]
          # create dataframe
          chrom_sizes <- data.frame(Chr,chrLengths[[input$chr]])
          colnames(chrom_sizes) <- c("chromosome", "size")
          chrseq = getSeq(BSGENOME_SELECTED, Chr,
                          start = Start, end = End)
          name2 <- paste(AT_Table_data$chr, "-", AT_Table_data$Label)
          meta2 <- data.frame(name2  , gieStain2 = AT_Table_data$Label)
          bands2 <- data.frame(chr = input$chr,
                               start = AT_Table_data$StartPos,
                               end = AT_Table_data$EndPos)
          all_regions = GRanges(bands2$chr, IRanges(bands2$start, bands2$end), name = meta2$name2, gieStain = meta2$gieStain2)
          tiles <- tile(x = all_regions, width = nuclFreq)
          tiles <- unlist(tiles)
        }
        if (input$refGenome=="hg38") {
          refSeq <- "hg38"
          BSGENOME_SELECTED = BSgenome.Hsapiens.UCSC.hg38
          # assign chromosome from selection
          chrLengths = seqlengths(BSGENOME_SELECTED)
          Chr = input$chr
          Start = 1
          End = chrLengths[[input$chr]]
          # create dataframe
          chrom_sizes <- data.frame(Chr,chrLengths[[input$chr]])
          colnames(chrom_sizes) <- c("chromosome", "size")
          chrseq = getSeq(BSGENOME_SELECTED, Chr,
                          start = Start, end = End)
          name2 <- paste(AT_Table_data$chr, "-", AT_Table_data$Label)
                      meta2 <- data.frame(name2  , gieStain2 = AT_Table_data$Label)
                      bands2 <- data.frame(chr = input$chr,
                                          start = AT_Table_data$StartPos,
                                          end = AT_Table_data$EndPos)
                      all_regions = GRanges(bands2$chr, IRanges(bands2$start, bands2$end), name = meta2$name2, gieStain = meta2$gieStain2)
                      tiles <- tile(x = all_regions, width = nuclFreq)
                      tiles <- unlist(tiles)
          
        }
        if (input$refGenome=="T2T") {
          refSeq <- "T2T"
          fasta_file = paste0(getwd(),"/Data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna")
          fa = FaFile(fasta_file) 
          BSGENOME_SELECTED = as(seqinfo(fa), "GRanges")
          # assign chromosome from selection
          chrLengths = BSGENOME_SELECTED@seqinfo@seqlengths
          Chr = input$chr
          Start = 1
          End = chrLengths[strtoi(str_remove(Chr, "chr"), base=0L)]
          # create dataframe
          chrom_sizes <- data.frame(Chr,End)
          colnames(chrom_sizes) <- c("chromosome", "size")
          Chr_name_t2t = BSGENOME_SELECTED@seqinfo@seqnames
          Chr_name_t2t = Chr_name_t2t[strtoi(str_remove(Chr, "chr"), base=0L)]
          # get chromosome sequence
          chrseq = getSeq(fa, BSGENOME_SELECTED[strtoi(str_remove(Chr, "chr"), base=0L)],
                          start = Start, end = End)
          tiles <- tileGenome(seqinfo(chrseq), tilewidth=nuclFreq)
          tiles <- unlist(tiles)
        }
        
        if (file.exists("AT_Table.txt")) {
          AT_Table_data = read.table('AT_Table.txt', header = TRUE, sep = "\t", dec = ".")
          print(AT_Table_data)
          AT_Table_data_subset = AT_Table_data[,1:4]
          colnames(AT_Table_data_subset) <- c("chromosome","start","end","Label")
          
          
          
          FS = filter(AT_Table_data_subset,Label!="control")
          p1<-ggplot()+
            geom_rect(data=chrom_sizes,
                      aes(xmin = 0, xmax = size, 
                          ymin = 0, ymax = 1),
                      colour="black", fill = "white")+
            ylim(0,10) +
            # black & white color theme 
            theme(axis.text.x = element_text(colour = "black"), 
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(), 
                  panel.background = element_blank(),
                  axis.ticks.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.title.y=element_blank()) +
            
            # add bands for FS
            geom_rect(data = FS, aes(xmin =start, 
                                              xmax = end, 
                                              ymax = 1, ymin = 0),
                      colour="black", fill = "grey")+
            geom_text_repel(data = FS, 
                            aes(y =1.2, x = start, label = Label, angle = 90,size=1), 
                            color = "black", show.legend = FALSE,nudge_y = 1.2)+ scale_radius(range = c(1,4))+
            ggtitle(paste(Chr,"Fragile Sites"))+
          # supress scientific notation on the y-axis
          scale_x_continuous(labels = comma) +
            xlab("region (bp)")
          
          # calculate Density per nuclFreq
          
          
          tilesDF = data.frame(tiles)
          
          RangeStart = tilesDF$start
          RangeEnd = tilesDF$end
          char_chrseq = toString(chrseq)
          
          
          countsDF = tilesDF
          countsDF$Freq_A = 0
          countsDF$Freq_T = 0
          countsDF$Freq_C = 0
          countsDF$Freq_G = 0
          withProgress(message = 'Processing', value = 0, {
          for (i in 1:nrow(countsDF)) {
            incProgress(1/i, detail = paste("Computing AT content per Interval ", i))
            
            s = RangeStart[i]
            e = RangeEnd[i]
            subsequ = substr(char_chrseq, s, e)
            w = tilesDF[i,4]
            
            Freq_A = lengths(regmatches(subsequ, gregexpr("A", subsequ)))
            Freq_T = lengths(regmatches(subsequ, gregexpr("T", subsequ)))
            Freq_C = lengths(regmatches(subsequ, gregexpr("C", subsequ)))
            Freq_G = lengths(regmatches(subsequ, gregexpr("G", subsequ)))
            
            countsDF$Freq_A[i] = Freq_A
            countsDF$Freq_T[i] = Freq_T
            countsDF$Freq_C[i] = Freq_C
            countsDF$Freq_G[i] = Freq_G
          }
          })
          # percent_GC = countsDF$Freq_G+countsDF$Freq_C
          # percent_AT = countsDF$Freq_A+countsDF$Freq_T
          percent_GC = (countsDF$Freq_G+countsDF$Freq_C)/(countsDF$Freq_A+countsDF$Freq_T+countsDF$Freq_C+countsDF$Freq_G)
          percent_AT = 1 - percent_GC
          countsDF$percent_GC = percent_GC
          countsDF$percent_AT = percent_AT
          p2<-ggplot() +
            geom_line(data=countsDF, aes(x=start, y=percent_AT),color="black",size=0.5) +
            scale_color_brewer(palette="Paired")+
            theme_minimal()+ylim(0,1)+
            ylab("GC Content")+scale_x_continuous(labels = comma)
          
          
          p2/p1
        }
 
      }
    })
    
    # output$plot_karyo <- renderPlot({
    #     
    #     # make sure input$nuclFreqInterval is between 1000 and 1000000
    #     nuclFreq = input_nuclFreqInterval()
    #     if (!between(nuclFreq, 1000, 1000000)) {
    #         showNotification("Nucleotide frequency interval not between 1000 and 1000000 bp!")
    #     }
    #     if (between(nuclFreq, 1000, 1000000)) {
    #         # calculate frequency of nucleotides per interval
    #       if (input$refGenome!="T2T") {
    #         kp <- plotKaryotype(plot.type = 2, chromosomes = input$chr,genome = input$refGenome)
    #       } else {
    #         refSeq == "T2T"
    #         fasta_file = paste0(getwd(),"/Data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna")
    #         fa = FaFile(fasta_file) 
    #         BSGENOME_SELECTED = as(seqinfo(fa), "GRanges")
    #         chrLengths = BSGENOME_SELECTED@seqinfo@seqlengths
    #         Chr = input$chr
    #         Start = 1
    #         End = chrLengths[strtoi(str_remove(Chr, "chr"), base=0L)]
    #         Chr_name_t2t = BSGENOME_SELECTED@seqinfo@seqnames
    #         Chr_name_t2t = Chr_name_t2t[strtoi(str_remove(Chr, "chr"), base=0L)]
    #         df_t2t = data.frame(chr=Chr,start=Start,end=End)
    #         df_t2t_gr = makeGRangesFromDataFrame(df_t2t) 
    #         # plot kp for t2t
    #         kp <- plotKaryotype(plot.type = 2, chromosomes = Chr,genome = df_t2t_gr)
    #       }
    #         
    #         # load regions from AT_Table if table exits
    #         #if (exists('df_AT_Table') && is.data.frame(get('df_AT_Table'))) {
    #             #AT_Table_data = df_AT_Table
    #         if (file.exists("AT_Table.txt")) {
    #             AT_Table_data = read.table('AT_Table.txt', header = TRUE, sep = "\t", dec = ".")
    #             # zoom in on detailed regions
    #             pos1 = min(AT_Table_data$StartPos)
    #             pos2 = max(AT_Table_data$EndPos)
    #             detail.region <- toGRanges(data.frame(input$chr, pos1, pos2))
    #             
    #             if (input$refGenome!="T2T") {
    #               kp <- plotKaryotype(plot.type = 2, zoom=detail.region,genome = input$refGenome)
    #             } else {
    #               fasta_file = paste0(getwd(),"/Data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna")
    #               fa = FaFile(fasta_file) 
    #               BSGENOME_SELECTED = as(seqinfo(fa), "GRanges")
    #               chrLengths = BSGENOME_SELECTED@seqinfo@seqlengths
    #               Chr = input$chr
    #               Start = 1
    #               End = chrLengths[strtoi(str_remove(Chr, "chr"), base=0L)]
    #               Chr_name_t2t = BSGENOME_SELECTED@seqinfo@seqnames
    #               Chr_name_t2t = Chr_name_t2t[strtoi(str_remove(Chr, "chr"), base=0L)]
    #               df_t2t = data.frame(chr=Chr,start=Start,end=End)
    #               df_t2t_gr = makeGRangesFromDataFrame(df_t2t) 
    #               
    #               pos1 = min(AT_Table_data$StartPos)
    #               pos2 = max(AT_Table_data$EndPos)
    #               detail.region <- toGRanges(data.frame(input$chr, pos1, pos2))
    #               
    #               kp <- plotKaryotype(plot.type = 2, zoom=detail.region,genome = df_t2t_gr)
    #             }
    #             
    #             #kpPlotGenes(kp, data = TxDb.Hsapiens.UCSC.hg19.knownGene, r0=0, r1=0.2, add.gene.names = FALSE, add.transcript.names = FALSE, add.strand.marks = FALSE)
    #             
    #             # create formal class GRanges
    #             nonControlRegions = AT_Table_data[AT_Table_data$Label!="control",]
    #             name <- paste(nonControlRegions$chr, "-", nonControlRegions$Label)
    #             meta <- data.frame(name  , gieStain = nonControlRegions$Label)
    #             bands <- data.frame(chr = input$chr, 
    #                                 start = nonControlRegions$StartPos, 
    #                                 end = nonControlRegions$EndPos)
    #             nonControl_regions <- GRanges(bands$chr, IRanges(bands$start, bands$end), name = meta$name, gieStain = meta$gieStain)
    #             # plot panel 1
    #             kpRect(kp, data = nonControl_regions, y0=0, y1=1, col="#FFDDDD", border=NA, r0=0, r1=0.8,clipping=TRUE,data.panel = 1)
    #             label <- nonControlRegions$Label
    #             x_pos = nonControlRegions$StartPos
    #             kpText(kp, chr=input$chr, x=x_pos, y=0.5, labels=label,srt=90,cex=0.8,data.panel = 1)
    #             kpAbline(kp,v = nonControlRegions$EndPos,r0=0, r1=0.8,lty=2,col="red")
    #             kpAbline(kp,v = nonControlRegions$StartPos,r0=0, r1=0.8,lty=2,col="red")
    #             # plot panel 2
    #             name2 <- paste(AT_Table_data$chr, "-", AT_Table_data$Label)
    #             meta2 <- data.frame(name2  , gieStain2 = AT_Table_data$Label)
    #             bands2 <- data.frame(chr = input$chr, 
    #                                 start = AT_Table_data$StartPos, 
    #                                 end = AT_Table_data$EndPos)
    #             all_regions = GRanges(bands2$chr, IRanges(bands2$start, bands2$end), name = meta2$name2, gieStain = meta2$gieStain2)
    #             tiles <- tile(x = all_regions, width = nuclFreq)
    #             tiles <- unlist(tiles)
    #             
    #             if (input$refGenome=="hg19") {
    #                 seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, tiles)
    #             }
    #             if (input$refGenome=="hg38") {
    #                 seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, tiles)
    #             }
    #             if (input$refGenome=="T2T") {
    #               refSeq == "T2T"
    #               fasta_file = paste0(getwd(),"/Data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna")
    #               fa = FaFile(fasta_file) 
    #               BSGENOME_SELECTED = as(seqinfo(fa), "GRanges")
    #               chrLengths = BSGENOME_SELECTED@seqinfo@seqlengths
    #               Chr = input$chr
    #               Start = 1
    #               End = chrLengths[strtoi(str_remove(Chr, "chr"), base=0L)]
    #               Chr_name_t2t = BSGENOME_SELECTED@seqinfo@seqnames
    #               Chr_name_t2t = Chr_name_t2t[strtoi(str_remove(Chr, "chr"), base=0L)]
    #               df_t2t = data.frame(chr=Chr,start=Start,end=End)
    #               df_t2t_gr = makeGRangesFromDataFrame(df_t2t) 
    #               # get chromosome sequence
    #               seqs = getSeq(fa, BSGENOME_SELECTED[strtoi(str_remove(Chr, "chr"), base=0L)],
    #                               start = Start, end = End)
    #               tiles <- tile(x = df_t2t_gr, width = nuclFreq)
    #               tiles <- unlist(tiles)
    #               SS = tiles@ranges@start
    #               SE = tiles@ranges@start+ tiles@ranges@width-1
    #               seqs = getSeq(fa, BSGENOME_SELECTED[strtoi(str_remove(Chr, "chr"), base=0L)],
    #                               start = SS, end = SE)
    #               
    #             }
    #             
    #             counts <- alphabetFrequency(seqs, baseOnly=TRUE)
    #             freqs <- counts/rowSums(counts)
    #             mcols(tiles) <- DataFrame(freqs[,bases])
    #             withProgress(message = 'Processing', value = 0, {
    #             cum.sums <- lapply(seq_len(length(tiles)), 
    #                                function(i) {
    #                                    incProgress(1/i, detail = paste("Calculating nucleotide frequency ", i))
    #                                    return(cumsum(as.numeric(data.frame(mcols(tiles))[i,c(1:4)])))
    #                                })
    #             })
    #             cum.sums <- do.call(rbind, cum.sums)
    #             cum.sums <- data.frame(cum.sums)
    #             names(cum.sums) <- bases
    #             
    #             kpLines(kp, data=tiles, y=cum.sums$G, lty=3, r0=2,r1=3,col="black",pch=".", cex=2)
    #             kpText(kp, data=tiles[length(tiles)], y=cum.sums$G[nrow(cum.sums)], r0=2, labels = "%GC", pos=4, clipping=FALSE)
    #         }
    #     }
    #    
    # })
    
    observeEvent(input$button_MotifCounts, {
        # to include motifs again
        motifL = strtoi(input$motifL)
        motifs = permutations(4,motifL, bases, repeats.allowed = TRUE)
        motifs = as.data.frame(motifs)
        nvars = colnames(motifs)
        ncols = length(nvars)
        pat = ""
        withProgress(message = 'Processing', value = 0, {
            for (i in 1:ncols) {
                incProgress(1/i, detail = paste("Generating Motifs", i))
                pat = paste(pat,motifs[,i])
            }
        })
        pat = as.data.frame(pat)
        pat = apply(pat,2,function(x)gsub('\\s+', '',x))
        colnames(pat) = c("Motifs")
        MotifsOnly_Table = pat
        
        if (file.exists("AT_Table.txt")) {
            DT_Table = read.table('AT_Table.txt', header = TRUE, sep = "\t", dec = ".")
            
            controlRegions_DT = DT_Table[DT_Table$Label=="control",]
            size_controlRegions = controlRegions_DT$EndPos-controlRegions_DT$StartPos +1
            L_control = sum(size_controlRegions)
            
            
            name <- paste(controlRegions_DT$chr, "-", controlRegions_DT$Label)
            meta <- data.frame(name  , gieStain = controlRegions_DT$Label)
            bands <- data.frame(chr = input$chr,
                                start = controlRegions_DT$StartPos,
                                end = controlRegions_DT$EndPos)
            controlRegions <- GRanges(bands$chr, IRanges(bands$start, bands$end), name = meta$name, gieStain = meta$gieStain)
            
            # do the counts for nonControl regions
            NoncontrolRegions_DT = DT_Table[DT_Table$Label!="control",]
            size_NoncontrolRegions = NoncontrolRegions_DT$EndPos-NoncontrolRegions_DT$StartPos +1
            #L_Noncontrol = sum(size_NoncontrolRegions)
            
            
            name2 <- paste(NoncontrolRegions_DT$chr, "-", NoncontrolRegions_DT$Label)
            meta2 <- data.frame(name2  , gieStain = NoncontrolRegions_DT$Label)
            bands2 <- data.frame(chr = input$chr,
                                 start = NoncontrolRegions_DT$StartPos,
                                 end = NoncontrolRegions_DT$EndPos)
            NoncontrolRegions <- GRanges(bands2$chr, IRanges(bands2$start, bands2$end), name = meta2$name2, gieStain = meta2$gieStain)
            
            
            if (input$refGenome=="hg19") {
                seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, controlRegions)
                seqs_NoncontrolRegions <- getSeq(BSgenome.Hsapiens.UCSC.hg19, NoncontrolRegions)
                
                sequence_controlRegions=DNAStringSet(as.character(seqs))
                sequence_NoncontrolRegions=DNAStringSet(as.character(seqs_NoncontrolRegions))
            }
            if (input$refGenome=="hg38") {
                seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, controlRegions)
                seqs_NoncontrolRegions <- getSeq(BSgenome.Hsapiens.UCSC.hg38, NoncontrolRegions)
                
                sequence_controlRegions=DNAStringSet(as.character(seqs))
                sequence_NoncontrolRegions=DNAStringSet(as.character(seqs_NoncontrolRegions))
            }
            if (input$refGenome=="T2T") {
              refSeq <- "T2T"
              fasta_file = paste0(getwd(),"/Data/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna")
              fa = FaFile(fasta_file) 
              BSGENOME_SELECTED = as(seqinfo(fa), "GRanges")
              # assign chromosome from selection
              chrLengths = BSGENOME_SELECTED@seqinfo@seqlengths
              Chr = input$chr
              Start = 1
              End = chrLengths[strtoi(str_remove(Chr, "chr"), base=0L)]
              # create dataframe
              chrom_sizes <- data.frame(Chr,End)
              colnames(chrom_sizes) <- c("chromosome", "size")
              Chr_name_t2t = BSGENOME_SELECTED@seqinfo@seqnames
              Chr_name_t2t = Chr_name_t2t[strtoi(str_remove(Chr, "chr"), base=0L)]
              # get chromosome sequence
              controlReg_startend = data.frame(controlRegions@ranges)
              NoncontrolRegions_startend = data.frame(NoncontrolRegions@ranges)
              
              chrseq = getSeq(fa, BSGENOME_SELECTED[strtoi(str_remove(Chr, "chr"), base=0L)],
                              start = Start, end = End)
              char_chrseq = toString(chrseq)
              
              DF_control_Seq = controlReg_startend
              DF_control_Seq$seq = ""
              
              DF_Noncontrol_Seq = NoncontrolRegions_startend
              DF_Noncontrol_Seq$seq = ""
              
              for (k in 1:nrow(controlReg_startend)) {
                tempRegion = controlReg_startend[k,]
                subsequ_k = substr(char_chrseq, tempRegion[1],tempRegion[2])
                DF_control_Seq$seq[k] = subsequ_k
              }
              
              for (k in 1:nrow(NoncontrolRegions_startend)) {
                tempRegion = NoncontrolRegions_startend[k,]
                subsequ_k = substr(char_chrseq, tempRegion[1],tempRegion[2])
                DF_Noncontrol_Seq$seq[k] = subsequ_k
              }
              
              sequence_controlRegions = DNAStringSet(x=DF_control_Seq$seq)
              sequence_NoncontrolRegions = DNAStringSet(x=DF_Noncontrol_Seq$seq)
              
            }
            
            count_current_pat_controlRegions = matrix(ncol = nrow(controlRegions_DT), nrow = 4^motifL)
            count_current_pat_NoncontrolRegions = matrix(ncol = nrow(NoncontrolRegions_DT), nrow = 4^motifL)
            
            countReps_current_pat_controlRegions = matrix(ncol = nrow(controlRegions_DT), nrow = 4^motifL)
            countReps_current_pat_NoncontrolRegions = matrix(ncol = nrow(NoncontrolRegions_DT), nrow = 4^motifL)
            
            # patterns
            motifs = permutations(4,motifL, bases, repeats.allowed = TRUE)
            motifs = as.data.frame(motifs)
            nvars = colnames(motifs)
            ncols = length(nvars)
            pat = ""
            
            withProgress(message = 'Processing', value = 0, {
                for (i in 1:ncols) {
                    incProgress(1/i, detail = paste("Generating Motifs", i))
                    pat = paste(pat,motifs[,i])
                }
            })
            pat = as.data.frame(pat)
            pat = apply(pat,2,function(x)gsub('\\s+', '',x))
            colnames(pat) = c("Motifs")
            
            if (input$motifOverlap == "Yes") {
              # counts Control
              withProgress(message = 'Processing Control Regions', value = 0, {
                for (i in 1:length(sequence_controlRegions)) {
                  incProgress(1/i, detail = paste("Generating Counts", i))
                  temp_Seq = as.character(sequence_controlRegions[[i]])
                  for (j in 1:4^motifL) {
                    current_pat = pat[j]
                    countpattern = str_count(temp_Seq, paste0("(?=",current_pat,")"))
                    count_current_pat_controlRegions[j,i] = countpattern
                    
                    # count repeats control
                    pos_pattern = data.frame(str_locate_all(temp_Seq, paste0("(?=",current_pat,")")))
                    d_pos_pattern = diff(pos_pattern$start)
                    countReps_current_pat_controlRegions[j,i] = length(which(d_pos_pattern<(motifL+1)))
                  }
                }
              })
              
              # counts non Control
              withProgress(message = 'Processing Non Control Regions', value = 0, {
                for (i in 1:length(sequence_NoncontrolRegions)) {
                  incProgress(1/i, detail = paste("Generating Counts", i))
                  temp_Seq = as.character(sequence_NoncontrolRegions[[i]])
                  for (j in 1:4^motifL) {
                    current_pat = pat[j]
                    countpattern = str_count(temp_Seq, paste0("(?=",current_pat,")"))
                    count_current_pat_NoncontrolRegions[j,i] = countpattern
                    
                    # count repeats non control
                    pos_pattern = data.frame(str_locate_all(temp_Seq, paste0("(?=",current_pat,")")))
                    d_pos_pattern = diff(pos_pattern$start)
                    countReps_current_pat_NoncontrolRegions[j,i] = length(which(d_pos_pattern<(motifL+1)))
                  }
                }
                
              })
            }
            
            if (input$motifOverlap == "No") {
              # counts Control
              withProgress(message = 'Processing Control Regions', value = 0, {
                for (i in 1:length(sequence_controlRegions)) {
                  incProgress(1/i, detail = paste("Generating Counts", i))
                  temp_Seq = as.character(sequence_controlRegions[[i]])
                  for (j in 1:4^motifL) {
                    current_pat = pat[j]
                    countpattern = str_count(temp_Seq, current_pat)
                    count_current_pat_controlRegions[j,i] = countpattern
                    
                    # count repeats control
                    pos_pattern = data.frame(str_locate_all(temp_Seq, current_pat))
                    d_pos_pattern = diff(pos_pattern$start)
                    countReps_current_pat_controlRegions[j,i] = length(which(d_pos_pattern<(motifL+1)))
                  }
                }
              })
              
              # counts non Control
              withProgress(message = 'Processing Non Control Regions', value = 0, {
                for (i in 1:length(sequence_NoncontrolRegions)) {
                  incProgress(1/i, detail = paste("Generating Counts", i))
                  temp_Seq = as.character(sequence_NoncontrolRegions[[i]])
                  for (j in 1:4^motifL) {
                    current_pat = pat[j]
                    countpattern = str_count(temp_Seq, current_pat)
                    count_current_pat_NoncontrolRegions[j,i] = countpattern
                    
                    # count repeats non control
                    pos_pattern = data.frame(str_locate_all(temp_Seq, current_pat))
                    d_pos_pattern = diff(pos_pattern$start)
                    countReps_current_pat_NoncontrolRegions[j,i] = length(which(d_pos_pattern<(motifL+1)))
                  }
                }
                
              })
            }
            
            # Count Occurrences 
            count_current_pat_controlRegions = rowSums(count_current_pat_controlRegions)
            norm_count_current_pat_controlRegions = count_current_pat_controlRegions/L_control
            count_current_pat_controlRegions = as.data.frame(count_current_pat_controlRegions)
            norm_count_current_pat_controlRegions = as.data.frame(norm_count_current_pat_controlRegions)
            
            L_CFS = NoncontrolRegions_DT$EndPos-NoncontrolRegions_DT$StartPos+1
            norm_count_current_pat_NoncontrolRegions = setNames(data.frame(matrix(ncol = length(L_CFS), nrow = 4^motifL)), paste0("norm_",NoncontrolRegions_DT$Label))
            for (l in 1:length(L_CFS)) {
              tempCFS_counts = count_current_pat_NoncontrolRegions[,l]
              temp_norm_CFS_counts = tempCFS_counts/L_CFS[l]
              norm_count_current_pat_NoncontrolRegions[,l] = temp_norm_CFS_counts
            }
            
            # Count Repeats
            countReps_current_pat_controlRegions = rowSums(countReps_current_pat_controlRegions)
            norm_countReps_current_pat_controlRegions = countReps_current_pat_controlRegions/L_control
            countReps_current_pat_controlRegions = as.data.frame(countReps_current_pat_controlRegions)
            norm_countReps_current_pat_controlRegions = as.data.frame(norm_countReps_current_pat_controlRegions)
            
            norm_countReps_current_pat_NoncontrolRegions = setNames(data.frame(matrix(ncol = length(L_CFS), nrow = 4^motifL)), paste0("norm_",NoncontrolRegions_DT$Label))
            for (l in 1:length(L_CFS)) {
              tempCFS_counts = countReps_current_pat_NoncontrolRegions[,l]
              temp_norm_CFS_counts = tempCFS_counts/L_CFS[l]
              norm_countReps_current_pat_NoncontrolRegions[,l] = temp_norm_CFS_counts
            }
            
            
            Motifs_Table = cbind(count_current_pat_controlRegions,norm_count_current_pat_controlRegions,count_current_pat_NoncontrolRegions,norm_count_current_pat_NoncontrolRegions,
                                 countReps_current_pat_controlRegions,norm_countReps_current_pat_controlRegions,
                                 countReps_current_pat_NoncontrolRegions,norm_countReps_current_pat_NoncontrolRegions)
            NoncontrolRegions_Labels = as.vector(NoncontrolRegions_DT$Label)
            
            
            colnames(Motifs_Table) = c("Counts_ControlRegion","norm_Counts_ControlRegion",NoncontrolRegions_Labels,paste0("norm_",NoncontrolRegions_DT$Label),
                                       "Count_Repeats_ControlRegion","norm_Count_Repeats_ControlRegion",paste0("Count_Repeats_",NoncontrolRegions_DT$Label),
                                       paste0("norm_Count_Repeats_",NoncontrolRegions_DT$Label))
            
            df <- Motifs_Table
            df = cbind(MotifsOnly_Table,df)
            DT_Table = df
            
        }
        else {
            DT_Table = MotifsOnly_Table
        }
        # output$MotifCounts_Table <- DT::renderDataTable(xtensions = 'Buttons', options = list(
        #   dom = 'Bfrtip',
        #   buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
        # ),{DT::datatable(DT_Table))}
        # 
        # output$MotifCounts_Table <- DT::renderDataTable(extensions = 'Buttons', options = list(
        #   dom = 'Bfrtip',
        #   buttons = list(
        #     list(extend = "csv", text = "Download Current Page", filename = "page",
        #          exportOptions = list(
        #            modifier = list(page = "current")
        #          )
        #     ),
        #     list(extend = "csv", text = "Download Full Results", filename = "data",
        #          exportOptions = list(
        #            modifier = list(page = "all")
        #          )
        #     )
        #   )
        #   
        # ),{
        #   DT_Table
        # })
        output$MotifCounts_Table <- DT::renderDT(server = FALSE, {
          DT::datatable(
            DT_Table,
            extensions = c("Buttons"),
            options = list(
              dom = 'Bfrtip',
              buttons = list(
                list(extend = "csv", text = "Download Current Page", filename = "page",
                     exportOptions = list(
                       modifier = list(page = "current")
                     )
                ),
                list(extend = "csv", text = "Download Full Results", filename = "data",
                     exportOptions = list(
                       modifier = list(page = "all")
                     )
                )
              )
            )
          )
        })
        
    })
    
    # # MotifPositions_Table 
    # observeEvent(input$button_MotifPositions, {
    #     # to include motifs again
    #     motifL = strtoi(input$motifL)
    #     motifs = permutations(4,motifL, bases, repeats.allowed = TRUE)
    #     motifs = as.data.frame(motifs)
    #     nvars = colnames(motifs)
    #     ncols = length(nvars)
    #     pat = ""
    #     withProgress(message = 'Processing', value = 0, {
    #         for (i in 1:ncols) {
    #             incProgress(1/i, detail = paste("Generating Motifs", i))
    #             pat = paste(pat,motifs[,i])
    #         }
    #     })
    #     pat = as.data.frame(pat)
    #     pat = apply(pat,2,function(x)gsub('\\s+', '',x))
    #     colnames(pat) = c("Motifs")
    #     MotifsOnly_Table = pat
    #     
    #     if (file.exists("AT_Table.txt")) {
    #         DT_Table = read.table('AT_Table.txt', header = TRUE, sep = "\t", dec = ".")
    #         
    #         controlRegions_DT = DT_Table[DT_Table$Label=="control",]
    #         size_controlRegions = controlRegions_DT$EndPos-controlRegions_DT$StartPos +1
    #         L_control = sum(size_controlRegions)
    #         
    #         
    #         name <- paste(controlRegions_DT$chr, "-", controlRegions_DT$Label)
    #         meta <- data.frame(name  , gieStain = controlRegions_DT$Label)
    #         bands <- data.frame(chr = input$chr,
    #                             start = controlRegions_DT$StartPos,
    #                             end = controlRegions_DT$EndPos)
    #         controlRegions <- GRanges(bands$chr, IRanges(bands$start, bands$end), name = meta$name, gieStain = meta$gieStain)
    #         
    #         # do the counts for nonControl regions
    #         NoncontrolRegions_DT = DT_Table[DT_Table$Label!="control",]
    #         size_NoncontrolRegions = NoncontrolRegions_DT$EndPos-NoncontrolRegions_DT$StartPos +1
    #         #L_Noncontrol = sum(size_NoncontrolRegions)
    #         
    #         
    #         name2 <- paste(NoncontrolRegions_DT$chr, "-", NoncontrolRegions_DT$Label)
    #         meta2 <- data.frame(name2  , gieStain = NoncontrolRegions_DT$Label)
    #         bands2 <- data.frame(chr = input$chr,
    #                              start = NoncontrolRegions_DT$StartPos,
    #                              end = NoncontrolRegions_DT$EndPos)
    #         NoncontrolRegions <- GRanges(bands2$chr, IRanges(bands2$start, bands2$end), name = meta2$name2, gieStain = meta2$gieStain)
    #         
    #         
    #         if (input$refGenome=="hg19") {
    #             seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, controlRegions)
    #             seqs_NoncontrolRegions <- getSeq(BSgenome.Hsapiens.UCSC.hg19, NoncontrolRegions)
    #         }
    #         if (input$refGenome=="hg38") {
    #             seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, controlRegions)
    #             seqs_NoncontrolRegions <- getSeq(BSgenome.Hsapiens.UCSC.hg38, NoncontrolRegions)
    #         }
    #         sequence_controlRegions=DNAStringSet(as.character(seqs))
    #         pos_current_pat_controlRegions = matrix(ncol = nrow(controlRegions_DT), nrow = 4^motifL)
    #         sequence_NoncontrolRegions=DNAStringSet(as.character(seqs_NoncontrolRegions))
    #         pos_current_pat_NoncontrolRegions = matrix(ncol = nrow(NoncontrolRegions_DT), nrow = 4^motifL)
    #         
    #         # patterns
    #         motifs = permutations(4,motifL, bases, repeats.allowed = TRUE)
    #         motifs = as.data.frame(motifs)
    #         nvars = colnames(motifs)
    #         ncols = length(nvars)
    #         pat = ""
    #         
    #         withProgress(message = 'Processing', value = 0, {
    #             for (i in 1:ncols) {
    #                 incProgress(1/i, detail = paste("Generating Motifs", i))
    #                 pat = paste(pat,motifs[,i])
    #             }
    #         })
    #         pat = as.data.frame(pat)
    #         pat = apply(pat,2,function(x)gsub('\\s+', '',x))
    #         colnames(pat) = c("Motifs")
    #         
    #         # positions Control
    #         withProgress(message = 'Processing Control Regions', value = 0, {
    #             for (i in 1:length(sequence_controlRegions)) {
    #                 incProgress(1/i, detail = paste("Generating Counts", i))
    #                 temp_Seq = as.character(sequence_controlRegions[[i]])
    #                 # temp Seq start and end pos
    #                 TEMP_START = controlRegions_DT$StartPos[i]
    #                 TEMP_END = controlRegions_DT$EndPos[i]
    #                 
    #                 for (j in 1:4^motifL) {
    #                     current_pat = pat[j]
    #                     Pos_pattern = str_locate_all(temp_Seq, current_pat)
    #                     Pos_pattern = as.data.frame(Pos_pattern)
    #                     Pos_pattern$start = Pos_pattern$start+TEMP_START-1
    #                     #print(Pos_pattern$start+TEMP_START-1)
    #                     #print(paste(Pos_pattern$start, collapse =","))
    #                     pos_current_pat_controlRegions[j,i] = paste(Pos_pattern$start, collapse =",")
    #                 }
    #             }
    #         })
    #         #print(pos_current_pat_controlRegions)
    #         # positions non Control
    #         withProgress(message = 'Processing Non Control Regions', value = 0, {
    #             for (i in 1:length(sequence_NoncontrolRegions)) {
    #                 incProgress(1/i, detail = paste("Generating Counts", i))
    #                 temp_Seq = as.character(sequence_NoncontrolRegions[[i]])
    #                 # temp Seq start and end pos
    #                 TEMP_START = NoncontrolRegions_DT$StartPos[i]
    #                 TEMP_END = NoncontrolRegions_DT$EndPos[i]
    #                 for (j in 1:4^motifL) {
    #                     current_pat = pat[j]
    #                     Pos_pattern = str_locate(temp_Seq, current_pat)
    #                     Pos_pattern = as.data.frame(Pos_pattern)
    #                     Pos_pattern$start = Pos_pattern$start+TEMP_START-1
    #                     pos_current_pat_NoncontrolRegions[j,i] = paste(Pos_pattern$start, collapse =",")
    #                 }
    #             }
    #         })
    #         
    #         Motifs_Table = cbind(pos_current_pat_controlRegions,pos_current_pat_NoncontrolRegions)
    #         NoncontrolRegions_Labels = as.vector(NoncontrolRegions_DT$Label)
    #         
    #         #colnames(Motifs_Table) = c("Positions_ControlRegion",NoncontrolRegions_Labels)
    #         df <- Motifs_Table
    #         df = cbind(MotifsOnly_Table,df)
    #         DT_Table = df
    #         
    #     }
    #     else{
    #         DT_Table = MotifsOnly_Table
    #     }
    #     
    #     output$MotifPositions_Table <- DT::renderDT({datatable(DT_Table,extensions = 'Buttons',
    #                                                                      options = list(scrollY = 300, scroller = TRUE, scrollX = T,
    #                                                                                              autoWidth = TRUE,
    #                                                                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
    #                                                                                     ))},server = FALSE)
    # })
    # 
    output$MotifCountsTableUpload <- DT::renderDataTable({
      userInputFile <- input$motifCountsTable
      ext <- tools::file_ext(userInputFile$datapath)
      req(userInputFile)
      
      # check file extenstion
      if (ext=="txt") {
        # read text file
        df <- read.table(userInputFile$datapath, header=TRUE, sep = "\t")
      }
      if (ext=="csv") {
        df <- read.table(userInputFile$datapath, header=TRUE, sep = ",")
      }
    })
    
    outVar = reactive({
      userInputFile <- input$motifCountsTable
      ext <- tools::file_ext(userInputFile$datapath)
      req(userInputFile)
      
      # check file extenstion
      if (ext=="txt") {
        # read text file
        df <- read.table(userInputFile$datapath, header=TRUE, sep = "\t")
      }
      if (ext=="csv") {
        df <- read.table(userInputFile$datapath, header=TRUE, sep = ",")
      }
      
      mydata = df
      colnames(df)
    })
    
    observe({
      updateSelectInput(session, "MotifNames",
                        choices = outVar()
      )})
    observe({
      updateSelectInput(session, "Counts_ControlRegion",
                        choices = outVar()
      )})
    observe({
      updateSelectInput(session, "norm_Counts_ControlRegion",
                        choices = outVar()
      )})
    observe({
      updateSelectInput(session, "Count_Repeats_ControlRegion",
                        choices = outVar()
      )})
    observe({
      updateSelectInput(session, "norm_Count_Repeats_ControlRegion",
                        choices = outVar()
      )})
    observe({
      updateSelectInput(session, "countROI",
                        choices = outVar()
      )})
    observe({
      updateSelectInput(session, "norm_countROI",
                        choices = outVar()
      )})
    observe({
      updateSelectInput(session, "repeat_countROI",
                        choices = outVar()
      )})
    observe({
      updateSelectInput(session, "norm_repeat_countROI",
                        choices = outVar()
      )})
    
    observeEvent(input$Stats_MotifsPresentHigherRate, {
      userInputFile <- input$motifCountsTable
      ext <- tools::file_ext(userInputFile$datapath)
      req(userInputFile)
      
      # check file extenstion
      if (ext=="txt") {
        # read text file
        df <- read.table(userInputFile$datapath, header=TRUE, sep = "\t")
      }
      if (ext=="csv") {
        df <- read.table(userInputFile$datapath, header=TRUE, sep = ",")
      }
      
      uploaded_counts_df = df

      # parameters needed: Counts_ControlRegion, norm_Counts_ControlRegion,
      # countROI, norm_countROI, length_Control, length_ROI

      # p = (Counts_nonCFS + Counts_CFS(:,i))./(L_nonCFS + length(fs_chr.StartPos(i):fs_chr.EndPos(i)));
      
      C1 = df[, input$Counts_ControlRegion]
      C2 = df[, input$countROI]

      L1 = input$sizeControl
      L2 = input$sizeROI
      
      normC1 = df[, input$norm_Counts_ControlRegion]
      normC2 = df[, input$norm_countROI]
      
      # p_total_occ = (C1+C2)/(L1+L2)
      # 
      # z_occ = (normC1-normC2)/sqrt(p_total_occ*(1-p_total_occ)*(1/L1 + 1/L2))
      # 
      # p_from_z_occ1 = pnorm(z_occ, 0, 1)
      # p_from_z_occ2 = (1-pnorm(z_occ, 0, 1))
      # p_from_z_occ = data.frame(p_from_z_occ1,p_from_z_occ2)
      # p_rom_z = 2*apply(p_from_z_occ, 1, FUN = min)
      # 
      # fdr_occ = p.adjust(p_rom_z, method = "BH")
      # print(fdr_occ)
      
      p_values_table = numeric()  
      fdr_values = numeric() 
      for (i in 1:nrow(df)) {
        c1 = C1[i]
        c2 = C2[i]
        
        normc1 = normC1[i]
        normc2 = normC2[i]
        
        # p = (c1+c2)/(L1+L2)
        # z = (normc1-normc2)/sqrt(p*(1-p*(1/L1 + 1/L2)))
        # p_from_z_1 = pnorm(z, 0, 1)
        # p_from_z_2 = 1-pnorm(z, 0, 1)
        # p_from_z = 2*min(p_from_z_1,p_from_z_2)
        
        res = prop.test(x = c(c1, c2), n = c(L1, L2), alternative = "less")
        p_values_table = c(p_values_table, res[["p.value"]])
        fdr = p.adjust(res[["p.value"]], method = "BH")
        
        fdr_values = c(fdr_values,fdr)
        
        
      }
      
    })
    
    observeEvent(input$Stats_MotifsRepeatedHigherRate, {
      userInputFile <- input$motifCountsTable
      ext <- tools::file_ext(userInputFile$datapath)
      req(userInputFile)
      
      # check file extenstion
      if (ext=="txt") {
        # read text file
        df <- read.table(userInputFile$datapath, header=TRUE, sep = "\t")
      }
      if (ext=="csv") {
        df <- read.table(userInputFile$datapath, header=TRUE, sep = ",")
      }
      
      uploaded_counts_df = df
      
      C1 = df[, input$Count_Repeats_ControlRegion]
      C2 = df[, input$repeat_countROI]
      
      L1 = input$sizeControl
      L2 = input$sizeROI
      
      normC1 = df[, input$norm_Count_Repeats_ControlRegion]
      normC2 = df[, input$norm_repeat_countROI]
      
      p_values_table = numeric() 
      fdr_values = numeric() 
      for (i in 1:nrow(df)) {
        c1 = C1[i]
        c2 = C2[i]
        normc1 = normC1[i]
        normc2 = normC2[i]
        
        # p = (c1+c2)/(L1+L2)
        # z = (normc1-normc2)/sqrt(p*(1-p)*(1/L1 + 1/L2))
        # 
        # p_from_z_1 = pnorm(z, 0, 1)
        # p_from_z_2 = 1-pnorm(z, 0, 1)
        # p_from_z = 2*min(p_from_z_1,p_from_z_2)
        # 
        # p_values_table = c(p_values_table,p_from_z)
        # fdr = p.adjust(p_from_z, method = "fdr")
        res = prop.test(x = c(c1, c2), n = c(L1, L2), alternative = "less", correct = TRUE)
        p_values_table = c(p_values_table, res[["p.value"]])
        fdr = p.adjust(res[["p.value"]], method = "fdr")
        fdr_values = c(fdr_values,fdr)
          
      }
      sig_FDR = which(fdr_values < 0.001)
      print(fdr_values)

      # p_values_table = numeric()  
      # fdr_values = numeric() 
      # for (i in 1:nrow(df)) {
      #   c1 = C1[i]
      #   c2 = C2[i]
      #   
      #   normc1 = normC1[i]
      #   normc2 = normC2[i]
      #   
      #   # res = prop.test(x = c(c1, c2), n = c(L1, L2), alternative = "two.sided")
      #   # p_values_table = c(p_values_table, res[["p.value"]])
      #   
      #   p = (c1+c2)/(L1+L2)
      #   z = (normc1-normc2)/sqrt(p*(1-p*(1/L1 + 1/L2)))
      #   p_from_z_1 = pnorm(z, 0, 1)
      #   p_from_z_2 = 1-pnorm(z, 0, 1)
      #   p_from_z = 2*min(p_from_z_1,p_from_z_2)
      #   p_values_table = c(p_values_table,p_from_z)
      #   
      #   fdr = p.adjust(p_from_z, method = "BH")
      #   
      #   fdr_values = c(fdr_values,fdr)
      #   
      # }
      # 
      # sig_FDR = which(fdr_values <= 0.01)
      # 
      # print(df[sig_FDR,])
      
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
