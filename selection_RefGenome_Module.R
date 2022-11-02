selection_RefGenome_ModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
  selectInput("chr", "Chromosome:",
              c("Chr1" = "chr1",
                "Chr2" = "chr2",
                "Chr3" = "chr3",
                "Chr4" = "chr4",
                "Chr5" = "chr5",
                "Chr6" = "chr6",
                "Chr7" = "chr7",
                "Chr8" = "chr8",
                "Chr9" = "chr9",
                "Chr10" = "chr10",
                "Chr11" = "chr11",
                "Chr12" = "chr12",
                "Chr13" = "chr13",
                "Chr14" = "chr14",
                "Chr15" = "chr15",
                "Chr16" = "chr16",
                "Chr17" = "chr17",
                "Chr18" = "chr18",
                "Chr19" = "chr19",
                "Chr20" = "chr20",
                "Chr21" = "chr21",
                "Chr22" = "chr22",
                "ChrX" = "chrX",
                "ChrY" = "chrY")),
  tableOutput('chrCoord')
  )
}

selection_RefGenome_Module <- function(input, output, session,chr_start_end) {
  # get chr info
  chr_coord <- reactive({
    filter(chr_start_end, chrom == input$chr)
  })
  
  output$chrCoord <- renderTable(chr_coord())
}