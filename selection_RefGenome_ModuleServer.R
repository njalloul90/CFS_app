selection_RefGenome_ModuleServer <- function(id) {
  moduleServer(id, 
               function(input, output, session) {
                 
                 selectedRef <- reactive({
                   input$refGenome
                 })
                 
                 return(selectedRef)
               }
               )
  
}