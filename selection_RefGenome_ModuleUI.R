selection_RefGenome_ModuleUI <- function(id) {
  ns <- NS(id)

  tagList(
  selectInput(ns("refGenome"), "Reference Genome:",
              c("GRCh37 hg19" = "hg19",
                "GRCh38 hg38" = "hg38")),
  textOutput(ns("selectedRef"))
  )
}

