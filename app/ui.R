#------------------------------------------------------------------------------
# PhyloCorrelations v1.0
# ui.R
# Last modified: 2020-02-22 12:33:18 (CET)
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
        "<a href='", CONFIGS$IndexURL, "'>PhyloCorrelate</a>"
      ))
    ),

    # profvis_ui("profiler"),

    tabPanel("KEGG Orthologs",
      useShinyjs(),
      div(
        id = "UI_KO_LOADING", h2("Loading..."),
        style = "z-index:1;position:absolute;background-color:white;width:100%;text-align:center;"
      ),
      uiOutput("UI_KO")
    ),

    tabPanel("TIGRFAM",
      useShinyjs(),
      div(
        id = "UI_TIGRFAM_LOADING", h2("Loading..."),
        style = "z-index:1;position:absolute;background-color:white;width:100%;text-align:center;"
      ),
      uiOutput("UI_TIGRFAM")
    ),

    tabPanel("PFAM",
      useShinyjs(),
      div(
        id = "UI_PFAM_LOADING", h2("Loading..."),
        style = "z-index:1;position:absolute;background-color:white;width:100%;text-align:center;"
      ),
      uiOutput("UI_PFAM")
    ),

    navbarMenu("BLASTP",

      tabPanel("Input & jobs",
        uiOutput("UI_BLASTP_INPUT")
      , value = "BLASTP_INPUT_TAB"),

      tabPanel("Results - KO",
        useShinyjs(),
        div(
          id = "UI_BLASTP_KO_LOADING", h2("Loading..."),
          style = "z-index:1;position:absolute;background-color:white;width:100%;text-align:center;"
        ),
        uiOutput("UI_BLASTP_KO")
      ),

      tabPanel("Results - TIGRFAM",
        useShinyjs(),
        div(
          id = "UI_BLASTP_TIGRFAM_LOADING", h2("Loading..."),
          style = "z-index:1;position:absolute;background-color:white;width:100%;text-align:center;"
        ),
        uiOutput("UI_BLASTP_TIGRFAM")
      ),

      tabPanel("Results - PFAM",
        useShinyjs(),
        div(
          id = "UI_BLASTP_PFAM_LOADING", h2("Loading..."),
          style = "z-index:1;position:absolute;background-color:white;width:100%;text-align:center;"
        ),
        uiOutput("UI_BLASTP_PFAM")
      )

    ),

    tabPanel("Help",
      uiOutput("UI_HELP")
    )

  )

)

msg("Finished loading ui.R")
ui
