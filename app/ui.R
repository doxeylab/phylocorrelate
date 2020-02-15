#------------------------------------------------------------------------------
# PhyloCorrelations v1.0
# ui.R
# Last modified: 2020-02-15 14:06:52 (CET)
# BJM Tremblay

msg("Loading ui.R")
ui <- function(request) fluidPage(

  theme = shinytheme("flatly"),

  title = "PhyloCorrelate: a database of phylogenetically correlating gene/protein families",

  br(),br(),br(),

  navbarPage(#"PhyloCorrelate",
    id = "NAVBARPG",
    position = "fixed-top",
    collapsible = TRUE,
    title = div(
      HTML(paste0(
        "<a href='", CONFIGS$URL, "'>PhyloCorrelate</a>"
      ))
    ),

    tabPanel("PFAM",
      uiOutput("UI_PFAM")
    ),

    tabPanel("TIGRFAM",
      uiOutput("UI_TIGRFAM")
    ),

    tabPanel("KEGG Orthologs",
      uiOutput("UI_KO")
    ),

    navbarMenu("BLASTP",

      tabPanel("Input & jobs",
        uiOutput("UI_BLASTP_INPUT")
      , value = "BLASTP_INPUT_TAB"),

      tabPanel("Results - PFAM",
        uiOutput("UI_BLASTP_PFAM")
      ),

      tabPanel("Results - TIGRFAM",
        uiOutput("UI_BLASTP_TIGRFAM")
      ),

      tabPanel("Results - KO",
        uiOutput("UI_BLASTP_KO")
      )
      
    ),

    tabPanel("Help",
      uiOutput("UI_HELP")
    )

  )

)

msg("Finished loading ui.R")
ui
