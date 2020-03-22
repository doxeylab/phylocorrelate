#------------------------------------------------------------------------------
# PhyloCorrelations v1.0
# server.R
# Last modified: 2020-03-22 22:43:48 (CET)
# BJM Tremblay

msg("Loading server.R")
server <- function(input, output, session) {

  # callModule(profvis_server, "profiler")
  # Rprof(
  #   strftime(Sys.time(), "%Y-%m-%d-%H-%M-%S.Rprof"),
  #   interval = 0.01, line.profiling = TRUE,
  #   gc.profiling = TRUE, memory.profiling = TRUE
  # )

  msg("Session start:", session$token)
  onSessionEnded({function() {
      msg("Session stop:", session$token)
      # Rprof(NULL)
    }
  })

  THRESH_MSG_ID <- character()

  #-----------------------------------------------------------------------------
  # Main Tab UI

  output$UI_PFAM <- renderUI(make_tab_main("PFAM"))
  output$UI_TIGRFAM <- renderUI(make_tab_main("TIGRFAM"))
  output$UI_KO <- renderUI(make_tab_main("KO"))
  output$UI_BLASTP_INPUT <- renderUI(make_tab_blastp_input())
  output$UI_BLASTP_PFAM <- renderUI(make_tab_main("PFAM", TRUE))
  output$UI_BLASTP_TIGRFAM <- renderUI(make_tab_main("TIGRFAM", TRUE))
  output$UI_BLASTP_KO <- renderUI(make_tab_main("KO",  TRUE))
  output$UI_HELP <- renderUI(make_tab_help())

  #-----------------------------------------------------------------------------
  # Side Input Panel

  SidePanelInputValues <- reactiveValues(
    INPUT_PFAM = "",
    INPUT_TIGRFAM = "",
    INPUT_KO = ""
  )

  output$SIDE_PANEL_INPUT_INFO_PFAM <- renderText({
    make_side_panel_input_info("PFAM",, SidePanelInputValues)
  })
  output$SIDE_PANEL_INPUT_INFO_TIGRFAM <- renderText({
    make_side_panel_input_info("TIGRFAM",, SidePanelInputValues)
  })
  output$SIDE_PANEL_INPUT_INFO_KO <- renderText({
    make_side_panel_input_info("KO",, SidePanelInputValues)
  })

  output$SIDE_PANEL_INPUT_SEARCH_TEXT_PFAM <- renderText({
    req(SidePanelInputValues$INPUT_PFAM)
    "Click on a PFAM to use it as input."
  })
  output$SIDE_PANEL_INPUT_SEARCH_TEXT_TIGRFAM <- renderText({
    req(SidePanelInputValues$INPUT_TIGRFAM)
    "Click on a TIGRFAM to use it as input."
  })
  output$SIDE_PANEL_INPUT_SEARCH_TEXT_KO <- renderText({
    req(SidePanelInputValues$INPUT_KO)
    "Click on a KO to use it as input."
  })

  output$SIDE_PANEL_INPUT_SEARCH_TABLE_PFAM <- renderDataTable({
    req(input$SIDE_PANEL_INPUT_SEARCH_PFAM)
    make_side_panel_table(searchPFAMs(input$SIDE_PANEL_INPUT_SEARCH_PFAM))
  })
  output$SIDE_PANEL_INPUT_SEARCH_TABLE_TIGRFAM <- renderDataTable({
    req(input$SIDE_PANEL_INPUT_SEARCH_TIGRFAM)
    make_side_panel_table(searchTIGRFAMs(input$SIDE_PANEL_INPUT_SEARCH_TIGRFAM))
  })
  output$SIDE_PANEL_INPUT_SEARCH_TABLE_KO <- renderDataTable({
    req(input$SIDE_PANEL_INPUT_SEARCH_KO)
    make_side_panel_table(searchKEGGs(input$SIDE_PANEL_INPUT_SEARCH_KO))
  })

  observe({
    lapply(c("_PFAM", "_TIGRFAM", "_KO"),
      function(x) {
        observeEvent(input[[p0(p0("SIDE_PANEL_INPUT_SEARCH_TABLE", x), "_cell_clicked")]], {
          info <- input[[p0(p0("SIDE_PANEL_INPUT_SEARCH_TABLE", x), "_cell_clicked")]]
          req(info$value)
          if (info$col != 0) return()
          SidePanelInputValues[[p0("INPUT", x)]] <- info$value
        })
      }
    )
  })

  output$SIDE_PANEL_INPUT_INFO_BLASTP_PFAM <- renderText({
    make_blastp_info_panel(BLASTP_LOADED$Id, BLASTP_LOADED$Info)
  })
  output$SIDE_PANEL_INPUT_INFO_BLASTP_TIGRFAM <- renderText({
    make_blastp_info_panel(BLASTP_LOADED$Id, BLASTP_LOADED$Info)
  })
  output$SIDE_PANEL_INPUT_INFO_BLASTP_KO <- renderText({
    make_blastp_info_panel(BLASTP_LOADED$Id, BLASTP_LOADED$Info)
  })

  #-----------------------------------------------------------------------------
  # Correlation Table

  observe({
    query <- parseQueryString(session$clientData$url_search)
    req(query[["goto"]])
    switch(substr(query[["goto"]], 1, 1),
      P = {
        if (query[["goto"]] %in% names(PFAMDesc)) {
          SidePanelInputValues$INPUT_PFAM <- query[["goto"]]
          updateNavbarPage(session, "NAVBARPG", "PFAM")
        } else {
          showModal(modalDialog(
            title = "Error",
            paste("Incorrect query,", query[["goto"]], "does not exist.")
          ))
        }
      },
      T = {
        if (query[["goto"]] %in% names(TIGRFAMDesc)) {
          SidePanelInputValues$INPUT_TIGRFAM <- query[["goto"]]
          updateNavbarPage(session, "NAVBARPG", "TIGRFAM")
        } else {
          showModal(modalDialog(
            title = "Error",
            paste("Incorrect query,", query[["goto"]], "does not exist.")
          ))
        }
      },
      K = {
        if (query[["goto"]] %in% names(KODesc)) {
          SidePanelInputValues$INPUT_KO <- query[["goto"]]
          updateNavbarPage(session, "NAVBARPG", "KEGG Orthologs")
        } else {
          showModal(modalDialog(
            title = "Error",
            paste("Incorrect query,", query[["goto"]], "does not exist.")
          ))
        }
      }
    )
  })

  output$MAIN_CORR_TABLE_PFAM <- DT::renderDataTable({
    req(SidePanelInputValues$INPUT_PFAM)
    out <- make_corr_table(SidePanelInputValues$INPUT_PFAM, "PFAM",
      make_list_globals("PFAM", input), THRESH_MSG_ID
    )
    THRESH_MSG_ID <- check_thresh_msg(out$THRESH_MSG_ID, out$keep)
    out$table
  })
  output$MAIN_CORR_TABLE_TIGRFAM <- DT::renderDataTable({
    req(SidePanelInputValues$INPUT_TIGRFAM)
    out <- make_corr_table(SidePanelInputValues$INPUT_TIGRFAM, "PFAM",
      make_list_globals("TIGRFAM", input), THRESH_MSG_ID
    )
    THRESH_MSG_ID <<- check_thresh_msg(out$THRESH_MSG_ID, out$keep)
    out$table
  })
  output$MAIN_CORR_TABLE_KO <- DT::renderDataTable({
    req(SidePanelInputValues$INPUT_KO)
    out <- make_corr_table(SidePanelInputValues$INPUT_KO, "PFAM",
      make_list_globals("KO", input), THRESH_MSG_ID
    )
    THRESH_MSG_ID <<- check_thresh_msg(out$THRESH_MSG_ID, out$keep)
    out$table
  })

  output$MAIN_CORR_TABLE_BLASTP_PFAM <- DT::renderDataTable({
    req(BLASTP_LOADED$Id)
    out <- make_corr_table(BLASTP_LOADED$Dat$PFAM, "PFAM",
      make_list_globals("PFAM", input, blastp = TRUE), THRESH_MSG_ID, isBlastp = TRUE
    )
    THRESH_MSG_ID <<- check_thresh_msg(out$THRESH_MSG_ID, out$keep)
    out$table
  })
  output$MAIN_CORR_TABLE_BLASTP_TIGRFAM <- DT::renderDataTable({
    req(BLASTP_LOADED$Id)
    out <- make_corr_table(BLASTP_LOADED$Dat$TIGRFAM, "TIGRFAM",
      make_list_globals("TIGRFAM", input, blastp = TRUE), THRESH_MSG_ID, isBlastp = TRUE
    )
    THRESH_MSG_ID <<- check_thresh_msg(out$THRESH_MSG_ID, out$keep)
    out$table
  })
  output$MAIN_CORR_TABLE_BLASTP_KO <- DT::renderDataTable({
    req(BLASTP_LOADED$Id)
    out <- make_corr_table(BLASTP_LOADED$Dat$KO, "KO",
      make_list_globals("KO", input, blastp = TRUE), THRESH_MSG_ID, isBlastp = TRUE
    )
    THRESH_MSG_ID <<- check_thresh_msg(out$THRESH_MSG_ID, out$keep)
    out$table
  })

  observeEvent(input$MAIN_CORR_TABLE_PFAM_cell_clicked, {
    info <- input$MAIN_CORR_TABLE_PFAM_cell_clicked
    req(info$value)
    if (info$col == 0) SidePanelInputValues$INPUT_PFAM <- info$value
  })
  observeEvent(input$MAIN_CORR_TABLE_TIGRFAM_cell_clicked, {
    info <- input$MAIN_CORR_TABLE_TIGRFAM_cell_clicked
    req(info$value)
    if (info$col == 0) SidePanelInputValues$INPUT_TIGRFAM <- info$value
  })
  observeEvent(input$MAIN_CORR_TABLE_KO_cell_clicked, {
    info <- input$MAIN_CORR_TABLE_KO_cell_clicked
    req(info$value)
    if (info$col == 0) SidePanelInputValues$INPUT_KO <- info$value
  })

  #-----------------------------------------------------------------------------
  # Correlation Network

  output$MAIN_CORR_NETWORK_PFAM <- renderVisNetwork({
    req(SidePanelInputValues$INPUT_PFAM)
    globals <- make_list_globals("PFAM", input,
      SidePanelInputValues, TabFilterCorrNetworkValues
    )
    future(makeCorrNetwork(globals$entry, globals)) %...>% return()
  })
  output$MAIN_CORR_NETWORK_TIGRFAM <- renderVisNetwork({
    req(SidePanelInputValues$INPUT_TIGRFAM)
    globals <- make_list_globals("TIGRFAM", input,
      SidePanelInputValues, TabFilterCorrNetworkValues
    )
    future(makeCorrNetwork(globals$entry, globals)) %...>% return()
  })
  output$MAIN_CORR_NETWORK_KO <- renderVisNetwork({
    req(SidePanelInputValues$INPUT_KO)
    globals <- make_list_globals("KO", input,
      SidePanelInputValues, TabFilterCorrNetworkValues
    )
    future(makeCorrNetwork(globals$entry, globals)) %...>% return()
  })

  output$MAIN_CORR_NETWORK_BLASTP_PFAM <- renderVisNetwork({
    req(BLASTP_LOADED$Id)
    globals <- make_list_globals("PFAM", input,
      TabFilterCorrNetworkValues = TabFilterCorrNetworkValues,
      BLASTP_LOADED = BLASTP_LOADED, blastp = TRUE
    )
    future(makeCorrNetwork(globals$PFAM, globals, isBlastp = TRUE)) %...>% return()
  })
  output$MAIN_CORR_NETWORK_BLASTP_TIGRFAM <- renderVisNetwork({
    req(BLASTP_LOADED$Id)
    globals <- make_list_globals("TIGRFAM", input,
      TabFilterCorrNetworkValues = TabFilterCorrNetworkValues,
      BLASTP_LOADED = BLASTP_LOADED, blastp = TRUE
    )
    future(makeCorrNetwork(globals$TIGRFAM, globals, isBlastp = TRUE)) %...>% return()
  })
  output$MAIN_CORR_NETWORK_BLASTP_KO <- renderVisNetwork({
    req(BLASTP_LOADED$Id)
    globals <- make_list_globals("KO", input,
      TabFilterCorrNetworkValues = TabFilterCorrNetworkValues,
      BLASTP_LOADED = BLASTP_LOADED, blastp = TRUE
    )
    future(makeCorrNetwork(globals$KO, globals, isBlastp = TRUE)) %...>% return()
  })

  #-----------------------------------------------------------------------------
  # GO/Pathway Enrichment

  output$MAIN_GOPATH_ENRICH_PFAM <- DT::renderDataTable({
    req(SidePanelInputValues$INPUT_PFAM)
    globals <- make_list_globals("PFAM", input,
      SidePanelInputValues, TabFilterCorrNetworkValues,
      TabFilterGopathEnrichValues
    )
    future(make_main_go_enrich(globals, "PFAM")) %...>% return()
  })
  output$MAIN_GOPATH_ENRICH_TIGRFAM <- DT::renderDataTable({
    req(SidePanelInputValues$INPUT_TIGRFAM)
    globals <- make_list_globals("TIGRFAM", input,
      SidePanelInputValues, TabFilterCorrNetworkValues,
      TabFilterGopathEnrichValues
    )
    future(make_main_go_enrich(globals, "TIGRFAM")) %...>% return()
  })
  output$MAIN_GOPATH_ENRICH_KO <- DT::renderDataTable({
    req(SidePanelInputValues$INPUT_KO)
    globals <- make_list_globals("KO", input,
      SidePanelInputValues, TabFilterCorrNetworkValues,
      TabFilterGopathEnrichValues
    )
    future(make_main_go_enrich(globals, "KO")) %...>% return()
  })

  output$MAIN_GOPATH_ENRICH_BLASTP_PFAM <- DT::renderDataTable({
    req(BLASTP_LOADED$Id)
    globals <- make_list_globals("PFAM", input,
      TabFilterGopathEnrichValues = TabFilterGopathEnrichValues,
      BLASTP_LOADED = BLASTP_LOADED, blastp = TRUE
    )
    future(make_main_go_enrich(globals, "PFAM", isBlastp = TRUE)) %...>% return()
  })
  output$MAIN_GOPATH_ENRICH_BLASTP_TIGRFAM <- DT::renderDataTable({
    req(BLASTP_LOADED$Id)
    globals <- make_list_globals("TIGRFAM", input,
      TabFilterGopathEnrichValues = TabFilterGopathEnrichValues,
      BLASTP_LOADED = BLASTP_LOADED, blastp = TRUE
    )
    future(make_main_go_enrich(globals, "TIGRFAM", isBlastp = TRUE)) %...>% return()
  })
  output$MAIN_GOPATH_ENRICH_BLASTP_KO <- DT::renderDataTable({
    req(BLASTP_LOADED$Id)
    globals <- make_list_globals("KO", input,
      TabFilterGopathEnrichValues = TabFilterGopathEnrichValues,
      BLASTP_LOADED = BLASTP_LOADED, blastp = TRUE
    )
    future(make_main_go_enrich(globals, "KO", isBlastp = TRUE)) %...>% return()
  })

  #-----------------------------------------------------------------------------
  # Score Distribution

  output$MAIN_SCORE_DIST_PFAM <- renderPlotly({
    req(SidePanelInputValues$INPUT_PFAM)
    globals <- make_list_globals("PFAM",
      SidePanelInputValues = SidePanelInputValues,
      TabFilterScoreDistValues = TabFilterScoreDistValues
    )
    future(makeDistPlot(globals$entry, globals$type)) %...>% return()
  })
  output$MAIN_SCORE_DIST_TIGRFAM <- renderPlotly({
    req(SidePanelInputValues$INPUT_TIGRFAM)
    globals <- make_list_globals("TIGRFAM",
      SidePanelInputValues = SidePanelInputValues,
      TabFilterScoreDistValues = TabFilterScoreDistValues
    )
    future(makeDistPlot(globals$entry, globals$type)) %...>% return()
  })
  output$MAIN_SCORE_DIST_KO <- renderPlotly({
    req(SidePanelInputValues$INPUT_KO)
    globals <- make_list_globals("KO",
      SidePanelInputValues = SidePanelInputValues,
      TabFilterScoreDistValues = TabFilterScoreDistValues
    )
    future(makeDistPlot(globals$entry, globals$type)) %...>% return()
  })

  output$MAIN_SCORE_DIST_BLASTP_PFAM <- renderPlotly({
    req(BLASTP_LOADED$Id)
    globals <- make_list_globals("PFAM",
      TabFilterScoreDistValues = TabFilterScoreDistValues,
      BLASTP_LOADED = BLASTP_LOADED, blastp = TRUE
    )
    future(makeDistPlot(globals$PFAM, globals$type, isBlastp = TRUE)) %...>% return()
  })
  output$MAIN_SCORE_DIST_BLASTP_TIGRFAM <- renderPlotly({
    req(BLASTP_LOADED$Id)
    globals <- make_list_globals("TIGRFAM",
      TabFilterScoreDistValues = TabFilterScoreDistValues,
      BLASTP_LOADED = BLASTP_LOADED, blastp = TRUE
    )
    future(makeDistPlot(globals$TIGRFAM, globals$type, isBlastp = TRUE)) %...>% return()
  })
  output$MAIN_SCORE_DIST_BLASTP_KO <- renderPlotly({
    req(BLASTP_LOADED$Id)
    globals <- make_list_globals("KO",
      TabFilterScoreDistValues = TabFilterScoreDistValues,
      BLASTP_LOADED = BLASTP_LOADED, blastp = TRUE
    )
    future(makeDistPlot(globals$KO, globals$type, isBlastp = TRUE)) %...>% return()
  })

  #-----------------------------------------------------------------------------
  # Matrix Comparison

  observeEvent(input$PFAM_MAT_BUTTON, {
    req(input$PFAM_MAT_BUTTON)
    updateTabsetPanel(session, "MAIN_TABS_PFAM",
      selected = "TAB_MAIN_MATRIX_COMP_PFAM"
    )
    TabFilterMatrixCompValues$FILTER_MATRIX_COMP_COMPARISON_VALUE_PFAM <- gsub(
      "button_", "", input$PFAM_MAT_BUTTON, fixed = TRUE
    )
  })
  observeEvent(input$TIGRFAM_MAT_BUTTON, {
    req(input$TIGRFAM_MAT_BUTTON)
    updateTabsetPanel(session, "MAIN_TABS_TIGRFAM",
      selected = "TAB_MAIN_MATRIX_COMP_TIGRFAM"
    )
    TabFilterMatrixCompValues$FILTER_MATRIX_COMP_COMPARISON_VALUE_TIGRFAM <- gsub(
      "button_", "", input$TIGRFAM_MAT_BUTTON, fixed = TRUE
    )
  })
  observeEvent(input$KO_MAT_BUTTON, {
    req(input$KO_MAT_BUTTON)
    updateTabsetPanel(session, "MAIN_TABS_KO",
      selected = "TAB_MAIN_MATRIX_COMP_KO"
    )
    TabFilterMatrixCompValues$FILTER_MATRIX_COMP_COMPARISON_VALUE_KO <- gsub(
      "button_", "", input$KO_MAT_BUTTON, fixed = TRUE
    )
  })

  observeEvent(input$B_PFAM_MAT_BUTTON, {
    req(input$B_PFAM_MAT_BUTTON)
    updateTabsetPanel(session, "MAIN_TABS_BLASTP_PFAM",
      selected = "TAB_MAIN_MATRIX_COMP_BLASTP_PFAM"
    )
    TabFilterMatrixCompValues$FILTER_MATRIX_COMP_COMPARISON_VALUE_BLASTP_PFAM <- gsub(
      "button_", "", input$B_PFAM_MAT_BUTTON, fixed = TRUE
    )
  })
  observeEvent(input$B_TIGRFAM_MAT_BUTTON, {
    req(input$B_TIGRFAM_MAT_BUTTON)
    updateTabsetPanel(session, "MAIN_TABS_BLASTP_TIGRFAM",
      selected = "TAB_MAIN_MATRIX_COMP_BLASTP_TIGRFAM"
    )
    TabFilterMatrixCompValues$FILTER_MATRIX_COMP_COMPARISON_VALUE_BLASTP_TIGRFAM <- gsub(
      "button_", "", input$B_TIGRFAM_MAT_BUTTON, fixed = TRUE
    )
  })
  observeEvent(input$B_KO_MAT_BUTTON, {
    req(input$B_KO_MAT_BUTTON)
    updateTabsetPanel(session, "MAIN_TABS_BLASTP_KO",
      selected = "TAB_MAIN_MATRIX_COMP_BLASTP_KO"
    )
    TabFilterMatrixCompValues$FILTER_MATRIX_COMP_COMPARISON_VALUE_BLASTP_KO <- gsub(
      "button_", "", input$B_KO_MAT_BUTTON, fixed = TRUE
    )
  })

  output$MAIN_MATRIX_COMP_PFAM <- renderPlot({
    req(SidePanelInputValues$INPUT_PFAM)
    req(TabFilterMatrixCompValues$FILTER_MATRIX_COMP_COMPARISON_VALUE_PFAM)
    globals <- make_list_globals("PFAM",
      SidePanelInputValues = SidePanelInputValues,
      TabFilterMatrixCompValues = TabFilterMatrixCompValues
    )
    future(plotMatrixComp(globals$entry, globals$comparison, globals$runs)) %...>% return()
  })
  output$MAIN_MATRIX_COMP_TIGRFAM <- renderPlot({
    req(SidePanelInputValues$INPUT_TIGRFAM)
    req(TabFilterMatrixCompValues$FILTER_MATRIX_COMP_COMPARISON_VALUE_TIGRFAM)
    globals <- make_list_globals("TIGRFAM",
      SidePanelInputValues = SidePanelInputValues,
      TabFilterMatrixCompValues = TabFilterMatrixCompValues
    )
    future(plotMatrixComp(globals$entry, globals$comparison, globals$runs)) %...>% return()
  })
  output$MAIN_MATRIX_COMP_KO <- renderPlot({
    req(SidePanelInputValues$INPUT_KO)
    req(TabFilterMatrixCompValues$FILTER_MATRIX_COMP_COMPARISON_VALUE_KO)
    globals <- make_list_globals("KO",
      SidePanelInputValues = SidePanelInputValues,
      TabFilterMatrixCompValues = TabFilterMatrixCompValues
    )
    future(plotMatrixComp(globals$entry, globals$comparison, globals$runs)) %...>% return()
  })

  output$MAIN_MATRIX_COMP_BLASTP_PFAM <- renderPlot({
    req(BLASTP_LOADED$Id)
    req(TabFilterMatrixCompValues$FILTER_MATRIX_COMP_COMPARISON_VALUE_BLASTP_PFAM)
    globals <- make_list_globals("PFAM",
      TabFilterMatrixCompValues = TabFilterMatrixCompValues,
      BLASTP_LOADED = BLASTP_LOADED, blastp = TRUE
    )
    future(plotMatrixComp(
      globals$Vec, globals$comparison, globals$runs, isBlastp = TRUE, globals$Id
    )) %...>% return()
  })
  output$MAIN_MATRIX_COMP_BLASTP_TIGRFAM <- renderPlot({
    req(BLASTP_LOADED$Id)
    req(TabFilterMatrixCompValues$FILTER_MATRIX_COMP_COMPARISON_VALUE_BLASTP_TIGRFAM)
    globals <- make_list_globals("TIGRFAM",
      TabFilterMatrixCompValues = TabFilterMatrixCompValues,
      BLASTP_LOADED = BLASTP_LOADED, blastp = TRUE
    )
    future(plotMatrixComp(
      globals$Vec, globals$comparison, globals$runs, isBlastp = TRUE, globals$Id
    )) %...>% return()
  })
  output$MAIN_MATRIX_COMP_BLASTP_KO <- renderPlot({
    req(BLASTP_LOADED$Id)
    req(TabFilterMatrixCompValues$FILTER_MATRIX_COMP_COMPARISON_VALUE_BLASTP_KO)
    globals <- make_list_globals("KO",
      TabFilterMatrixCompValues = TabFilterMatrixCompValues,
      BLASTP_LOADED = BLASTP_LOADED, blastp = TRUE
    )
    future(plotMatrixComp(
      globals$Vec, globals$comparison, globals$runs, isBlastp = TRUE, globals$Id
    )) %...>% return()
  })

  #-----------------------------------------------------------------------------
  # Blastp Input Page

  BLASTP_PARAMETERS <- reactiveValues(
    evalue = 0.00001,
    pident = 30,
    qlenp = 70,
    slenp = 70
  )

  BLASTP_LOADED <- reactiveValues(
    Id = NULL,
    Info = NULL,
    Dat = NULL
  )

  observeEvent(input$BLASTP_FILTER_EVAL, {
    req(input$BLASTP_FILTER_EVAL)
    BLASTP_PARAMETERS$evalue <- input$BLASTP_FILTER_EVAL
  })
  observeEvent(input$BLASTP_FILTER_PIDENT, {
    req(input$BLASTP_FILTER_PIDENT)
    BLASTP_PARAMETERS$pident <- input$BLASTP_FILTER_PIDENT
  })
  observeEvent(input$BLASTP_FILTER_QLENP, {
    req(input$BLASTP_FILTER_QLENP)
    BLASTP_PARAMETERS$qlenp <- input$BLASTP_FILTER_QLENP
  })
  observeEvent(input$BLASTP_FILTER_SLENP, {
    req(input$BLASTP_FILTER_SLENP)
    BLASTP_PARAMETERS$slenp <- input$BLASTP_FILTER_SLENP
  })
  observeEvent(input$BLASTP_INPUT_BUTTON, {
    req(input$BLASTP_INPUT_SEQUENCE)
    if (!CONFIGS$UseBlastp) {
      showModal(modalDialog( title = "Error",
        paste(
          "BLASTP functionality has been disabled. Enable it by editing",
          "the config file."
        )
      ))
      return()
    }
    startBlast()
    Id <- getBlastpId(
      input$BLASTP_INPUT_SEQUENCE,
      input$BLASTP_INPUT_EMAIL,
      BLASTP_PARAMETERS$evalue,
      BLASTP_PARAMETERS$pident,
      BLASTP_PARAMETERS$qlenp,
      BLASTP_PARAMETERS$slenp,
      session$token
    )
    if (is.null(Id$Id)) {
      showModal(modalDialog(title = "BLASTP job submission", Id$status))
    } else {
      showModal(modalDialog(
        title = "BLASTP job submission", paste(Id$status, "JOB ID", Id$Id))
      )
    }
  })

  observe({
    query <- parseQueryString(session$clientData$url_search)
    req(query[["blastp"]])
    Id <- as.integer(query[["blastp"]])
    startBlast()
    if (is.na(Id)) {
      showModal(modalDialog(
        title = "BLASTP Job ID checker",
        "The job ID must be a number."
      ))
    } else {
      res <- checkBlastpId(Id)
      if (res$load) {
        BLASTP_LOADED$Id <- Id
        BLASTP_LOADED$Info <- getBlastpInfo(Id)
        BLASTP_LOADED$Dat <- readRDS(paste0("blastp/app/", Id))
        modal_text <- paste(
          "Your job is complete and the results loaded. Open one of the",
          "Results tabs to view them."
        )
        if (BLASTP_LOADED$Dat$version < BlastpVer) {
          modal_text <- paste(
            modal_text,
            "Warning: this submission is from an older version of the analysis",
            "pipeline. Some data may be absent. Re-submit your query for updated",
            "results."
          )
        }
        updateTabsetPanel(session, "NAVBARPG", selected = "BLASTP_TABS_INPUT")
        showModal(modalDialog(
          title = "BLASTP Job ID checker", modal_text
        ))
      } else {
        showModal(modalDialog(title = "BLASTP Job ID checker", res$status))
      }
    }
  })

  observeEvent(input$BLASTP_CHECK_JOB_BUTTON, {
    req(input$BLASTP_CHECK_JOB_ID)
    startBlast()
    Id <- as.integer(input$BLASTP_CHECK_JOB_ID)
    if (is.na(Id)) {
      showModal(modalDialog(
        title = "BLASTP Job ID checker",
        "The job ID must be a number."
      ))
    } else {
      res <- checkBlastpId(Id)
      if (res$load) {
        BLASTP_LOADED$Id <- Id
        BLASTP_LOADED$Info <- getBlastpInfo(Id)
        BLASTP_LOADED$Dat <- readRDS(paste0("blastp/app/", Id))
        modal_text <- paste(
          "Your job is complete and the results loaded. Open one of the",
          "Results tabs to view them."
        )
        if (BLASTP_LOADED$Dat$version < BlastpVer) {
          modal_text <- paste(
            modal_text,
            "Warning: this submission is from an older version of the analysis",
            "pipeline. Some data may be absent. Re-submit your query for updated",
            "results."
          )
        }
        showModal(modalDialog(
          title = "BLASTP Job ID checker", modal_text
        ))
      } else {
        showModal(modalDialog(title = "BLASTP Job ID checker", res$status))
      }
    }
  })

  #-----------------------------------------------------------------------------
  # Tab-specific Filters UI

  TabFilterMatrixCompValues <- reactiveValues(
    FILTER_MATRIX_COMP_TYPE_VALUE_PFAM = "Naive",
    FILTER_MATRIX_COMP_TYPE_VALUE_TIGRFAM = "Naive",
    FILTER_MATRIX_COMP_TYPE_VALUE_KO = "Naive",
    FILTER_MATRIX_COMP_TYPE_VALUE_BLASTP_PFAM = "Naive",
    FILTER_MATRIX_COMP_TYPE_VALUE_BLASTP_TIGRFAM = "Naive",
    FILTER_MATRIX_COMP_TYPE_VALUE_BLASTP_KO = "Naive",
    FILTER_MATRIX_COMP_COMPARISON_VALUE_PFAM = "",
    FILTER_MATRIX_COMP_COMPARISON_VALUE_TIGRFAM = "",
    FILTER_MATRIX_COMP_COMPARISON_VALUE_KO = "",
    FILTER_MATRIX_COMP_COMPARISON_VALUE_BLASTP_PFAM = "",
    FILTER_MATRIX_COMP_COMPARISON_VALUE_BLASTP_TIGRFAM = "",
    FILTER_MATRIX_COMP_COMPARISON_VALUE_BLASTP_KO = ""
  )

  TabFilterScoreDistValues <- reactiveValues(
    FILTER_SCORE_DIST_METRIC_VALUE_PFAM = "rJC",
    FILTER_SCORE_DIST_METRIC_VALUE_TIGRFAM = "rJC",
    FILTER_SCORE_DIST_METRIC_VALUE_KO = "rJC",
    FILTER_SCORE_DIST_METRIC_VALUE_BLASTP_PFAM = "rJC",
    FILTER_SCORE_DIST_METRIC_VALUE_BLASTP_TIGRFAM = "rJC",
    FILTER_SCORE_DIST_METRIC_VALUE_BLASTP_KO = "rJC"
  )

  TabFilterGopathEnrichValues <- reactiveValues(
    FILTER_GOPATH_ENRICH_ONTOLOGY_VALUE_PFAM = "BP",
    FILTER_GOPATH_ENRICH_MAX_FDR_VALUE_PFAM = 0.05,
    FILTER_GOPATH_ENRICH_MAX_FP_VALUE_PFAM = 0.05,
    FILTER_GOPATH_ENRICH_MIN_TOTAL_VALUE_PFAM = 3,
    FILTER_GOPATH_ENRICH_MIN_PRESENT_VALUE_PFAM = 2,
    FILTER_GOPATH_ENRICH_MIN_EXPECTED_VALUE_PFAM = 0,
    FILTER_GOPATH_ENRICH_ONTOLOGY_VALUE_TIGRFAM = "BP",
    FILTER_GOPATH_ENRICH_MAX_FDR_VALUE_TIGRFAM = 0.05,
    FILTER_GOPATH_ENRICH_MAX_FP_VALUE_TIGRFAM = 0.05,
    FILTER_GOPATH_ENRICH_MIN_TOTAL_VALUE_TIGRFAM = 3,
    FILTER_GOPATH_ENRICH_MIN_PRESENT_VALUE_TIGRFAM = 2,
    FILTER_GOPATH_ENRICH_MIN_EXPECTED_VALUE_TIGRFAM = 0,
    FILTER_GOPATH_ENRICH_MAX_FDR_VALUE_KO = 0.05,
    FILTER_GOPATH_ENRICH_MAX_FP_VALUE_KO = 0.05,
    FILTER_GOPATH_ENRICH_MIN_TOTAL_VALUE_KO = 3,
    FILTER_GOPATH_ENRICH_MIN_PRESENT_VALUE_KO = 2,
    FILTER_GOPATH_ENRICH_MIN_EXPECTED_VALUE_KO = 0,
    FILTER_GOPATH_ENRICH_ONTOLOGY_VALUE_BLASTP_PFAM = "BP",
    FILTER_GOPATH_ENRICH_MAX_FDR_VALUE_BLASTP_PFAM = 0.05,
    FILTER_GOPATH_ENRICH_MAX_FP_VALUE_BLASTP_PFAM = 0.05,
    FILTER_GOPATH_ENRICH_MIN_TOTAL_VALUE_BLASTP_PFAM = 3,
    FILTER_GOPATH_ENRICH_MIN_PRESENT_VALUE_BLASTP_PFAM = 2,
    FILTER_GOPATH_ENRICH_MIN_EXPECTED_VALUE_BLASTP_PFAM = 0,
    FILTER_GOPATH_ENRICH_ONTOLOGY_VALUE_BLASTP_TIGRFAM = "BP",
    FILTER_GOPATH_ENRICH_MAX_FDR_VALUE_BLASTP_TIGRFAM = 0.05,
    FILTER_GOPATH_ENRICH_MAX_FP_VALUE_BLASTP_TIGRFAM = 0.05,
    FILTER_GOPATH_ENRICH_MIN_TOTAL_VALUE_BLASTP_TIGRFAM = 3,
    FILTER_GOPATH_ENRICH_MIN_PRESENT_VALUE_BLASTP_TIGRFAM = 2,
    FILTER_GOPATH_ENRICH_MIN_EXPECTED_VALUE_BLASTP_TIGRFAM = 0,
    FILTER_GOPATH_ENRICH_MAX_FDR_VALUE_BLASTP_KO = 0.05,
    FILTER_GOPATH_ENRICH_MAX_FP_VALUE_BLASTP_KO = 0.05,
    FILTER_GOPATH_ENRICH_MIN_TOTAL_VALUE_BLASTP_KO = 3,
    FILTER_GOPATH_ENRICH_MIN_PRESENT_VALUE_BLASTP_KO = 2,
    FILTER_GOPATH_ENRICH_MIN_EXPECTED_VALUE_BLASTP_KO = 0
  )

  TabFilterCorrNetworkValues <- reactiveValues(
    FILTER_CORR_NETWORK_METRIC_VALUE_PFAM = "CS",
    FILTER_CORR_NETWORK_MIN_SECONDARY_EDGE_VALUE_PFAM = 0.5,
    FILTER_CORR_NETWORK_CS_SECONDARY_EDGE_PFAM = c("Low", "High", "Very high"),
    FILTER_CORR_NETWORK_METRIC_VALUE_TIGRFAM = "CS",
    FILTER_CORR_NETWORK_MIN_SECONDARY_EDGE_VALUE_TIGRFAM = 0.5,
    FILTER_CORR_NETWORK_CS_SECONDARY_EDGE_TIGRFAM = c("Low", "High", "Very high"),
    FILTER_CORR_NETWORK_METRIC_VALUE_KO = "CS",
    FILTER_CORR_NETWORK_MIN_SECONDARY_EDGE_VALUE_KO = 0.5,
    FILTER_CORR_NETWORK_CS_SECONDARY_EDGE_KO = c("Low", "High", "Very high"),
    FILTER_CORR_NETWORK_METRIC_VALUE_BLASTP_PFAM = "CS",
    FILTER_CORR_NETWORK_MIN_SECONDARY_EDGE_VALUE_BLASTP_PFAM = 0.5,
    FILTER_CORR_NETWORK_CS_SECONDARY_EDGE_BLASTP_PFAM = c("Low", "High", "Very high"),
    FILTER_CORR_NETWORK_METRIC_VALUE_BLASTP_TIGRFAM = "CS",
    FILTER_CORR_NETWORK_MIN_SECONDARY_EDGE_VALUE_BLASTP_TIGRFAM = 0.5,
    FILTER_CORR_NETWORK_CS_SECONDARY_EDGE_BLASTP_TIGRFAM = c("Low", "High", "Very high"),
    FILTER_CORR_NETWORK_METRIC_VALUE_BLASTP_KO = "CS",
    FILTER_CORR_NETWORK_MIN_SECONDARY_EDGE_VALUE_BLASTP_KO = 0.5,
    FILTER_CORR_NETWORK_CS_SECONDARY_EDGE_BLASTP_KO = c("Low", "High", "Very high")
  )

  output$UI_SIDE_TAB_FILTERS_PFAM <- renderUI({
    switch(input$MAIN_TABS_PFAM,
      TAB_MAIN_CORR_TABLE_PFAM = make_tab_filter_empty(),
      TAB_MAIN_CORR_NETWORK_PFAM = make_tab_filter_network(
        "PFAM",, TabFilterCorrNetworkValues
      ),
      TAB_MAIN_GOPATH_ENRICH_PFAM = make_tab_filter_gopath(
        "PFAM",, TabFilterGopathEnrichValues
      ),
      TAB_MAIN_SCORE_DIST_PFAM = make_tab_filter_score_dist(
        "PFAM",, TabFilterScoreDistValues
      ),
      TAB_MAIN_MATRIX_COMP_PFAM = make_tab_filter_matrix_comp(
        "PFAM",, TabFilterMatrixCompValues
      )
    )
  })

  output$UI_SIDE_TAB_FILTERS_TIGRFAM <- renderUI({
    switch(input$MAIN_TABS_TIGRFAM,
      TAB_MAIN_CORR_TABLE_TIGRFAM = make_tab_filter_empty(),
      TAB_MAIN_CORR_NETWORK_TIGRFAM = make_tab_filter_network(
        "TIGRFAM",, TabFilterCorrNetworkValues
      ),
      TAB_MAIN_GOPATH_ENRICH_TIGRFAM = make_tab_filter_gopath(
        "TIGRFAM",, TabFilterGopathEnrichValues
      ),
      TAB_MAIN_SCORE_DIST_TIGRFAM = make_tab_filter_score_dist(
        "TIGRFAM",, TabFilterScoreDistValues
      ),
      TAB_MAIN_MATRIX_COMP_TIGRFAM = make_tab_filter_matrix_comp(
        "TIGRFAM",, TabFilterMatrixCompValues
      )
    )
  })

  output$UI_SIDE_TAB_FILTERS_KO <- renderUI({
    switch(input$MAIN_TABS_KO,
      TAB_MAIN_CORR_TABLE_KO = make_tab_filter_empty(),
      TAB_MAIN_CORR_NETWORK_KO = make_tab_filter_network(
        "KO",, TabFilterCorrNetworkValues
      ),
      TAB_MAIN_GOPATH_ENRICH_KO = make_tab_filter_gopath(
        "KO",, TabFilterGopathEnrichValues
      ),
      TAB_MAIN_SCORE_DIST_KO = make_tab_filter_score_dist(
        "KO",, TabFilterScoreDistValues
      ),
      TAB_MAIN_MATRIX_COMP_KO = make_tab_filter_matrix_comp(
        "KO",, TabFilterMatrixCompValues
      )
    )
  })

  output$UI_SIDE_TAB_FILTERS_BLASTP_PFAM <- renderUI({
    switch(input$MAIN_TABS_BLASTP_PFAM,
      TAB_MAIN_CORR_TABLE_BLASTP_PFAM = make_tab_filter_empty(),
      TAB_MAIN_CORR_NETWORK_BLASTP_PFAM = make_tab_filter_network(
        "PFAM", TRUE, TabFilterCorrNetworkValues
      ),
      TAB_MAIN_GOPATH_ENRICH_BLASTP_PFAM = make_tab_filter_gopath(
        "PFAM", TRUE, TabFilterGopathEnrichValues
      ),
      TAB_MAIN_SCORE_DIST_BLASTP_PFAM = make_tab_filter_score_dist(
        "PFAM", TRUE, TabFilterScoreDistValues
      ),
      TAB_MAIN_MATRIX_COMP_BLASTP_PFAM = make_tab_filter_matrix_comp(
        "PFAM", TRUE, TabFilterMatrixCompValues
      )
    )
  })

  output$UI_SIDE_TAB_FILTERS_BLASTP_TIGRFAM <- renderUI({
    switch(input$MAIN_TABS_BLASTP_TIGRFAM,
      TAB_MAIN_CORR_TABLE_BLASTP_TIGRFAM = make_tab_filter_empty(),
      TAB_MAIN_CORR_NETWORK_BLASTP_TIGRFAM = make_tab_filter_network(
        "TIGRFAM", TRUE, TabFilterCorrNetworkValues
      ),
      TAB_MAIN_GOPATH_ENRICH_BLASTP_TIGRFAM = make_tab_filter_gopath(
        "TIGRFAM", TRUE, TabFilterGopathEnrichValues
      ),
      TAB_MAIN_SCORE_DIST_BLASTP_TIGRFAM = make_tab_filter_score_dist(
        "TIGRFAM", TRUE, TabFilterScoreDistValues
      ),
      TAB_MAIN_MATRIX_COMP_BLASTP_TIGRFAM = make_tab_filter_matrix_comp(
        "TIGRFAM", TRUE, TabFilterMatrixCompValues
      )
    )
  })

  output$UI_SIDE_TAB_FILTERS_BLASTP_KO <- renderUI({
    switch(input$MAIN_TABS_BLASTP_KO,
      TAB_MAIN_CORR_TABLE_BLASTP_KO = make_tab_filter_empty(),
      TAB_MAIN_CORR_NETWORK_BLASTP_KO = make_tab_filter_network(
        "KO", TRUE, TabFilterCorrNetworkValues
      ),
      TAB_MAIN_GOPATH_ENRICH_BLASTP_KO = make_tab_filter_gopath(
        "KO", TRUE, TabFilterGopathEnrichValues
      ),
      TAB_MAIN_SCORE_DIST_BLASTP_KO = make_tab_filter_score_dist(
        "KO", TRUE, TabFilterScoreDistValues
      ),
      TAB_MAIN_MATRIX_COMP_BLASTP_KO = make_tab_filter_matrix_comp(
        "KO", TRUE, TabFilterMatrixCompValues
      )
    )
  })

  observe({
    lapply(
      c("_PFAM", "_TIGRFAM", "_KO", "_BLASTP_PFAM", "_BLASTP_TIGRFAM", "_BLASTP_KO"),
      function(x) {
        observeEvent(input[[p0("FILTER_MATRIX_COMP_TYPE", x)]], {
          req(input[[p0("FILTER_MATRIX_COMP_TYPE", x)]])
          TabFilterMatrixCompValues[[p0("FILTER_MATRIX_COMP_TYPE_VALUE", x)]] <- 
            input[[p0("FILTER_MATRIX_COMP_TYPE", x)]]
        })
        observeEvent(input[[p0("FILTER_SCORE_DIST_METRIC", x)]], {
          req(input[[p0("FILTER_SCORE_DIST_METRIC", x)]])
          TabFilterScoreDistValues[[p0("FILTER_SCORE_DIST_METRIC_VALUE", x)]] <- 
            input[[p0("FILTER_SCORE_DIST_METRIC", x)]]
        })
        observeEvent(input[[p0("FILTER_CORR_NETWORK_MIN_SECONDARY_EDGE", x)]], {
          req(input[[p0("FILTER_CORR_NETWORK_MIN_SECONDARY_EDGE", x)]])
          TabFilterCorrNetworkValues[[p0("FILTER_CORR_NETWORK_MIN_SECONDARY_EDGE_VALUE", x)]] <- 
            as.numeric(input[[p0("FILTER_CORR_NETWORK_MIN_SECONDARY_EDGE", x)]])
        })
        observeEvent(input[[p0("FILTER_CORR_NETWORK_METRIC", x)]], {
          req(input[[p0("FILTER_CORR_NETWORK_METRIC", x)]])
          TabFilterCorrNetworkValues[[p0("FILTER_CORR_NETWORK_METRIC_VALUE", x)]] <- 
            input[[p0("FILTER_CORR_NETWORK_METRIC", x)]]
        })
        observeEvent(input[[p0("FILTER_CORR_NETWORK_CS_SECONDARY_EDGE", x)]], {
          req(input[[p0("FILTER_CORR_NETWORK_CS_SECONDARY_EDGE", x)]])
          TabFilterCorrNetworkValues[[p0("FILTER_CORR_NETWORK_CS_SECONDARY_EDGE", x)]] <- input[[p0("FILTER_CORR_NETWORK_CS_SECONDARY_EDGE", x)]]
        })
        observeEvent(input[[p0("FILTER_GOPATH_ENRICH_MIN_EXPECTED", x)]], {
          req(input[[p0("FILTER_GOPATH_ENRICH_MIN_EXPECTED", x)]])
          TabFilterGopathEnrichValues[[p0("FILTER_GOPATH_ENRICH_MIN_EXPECTED_VALUE", x)]] <- 
            as.numeric(input[[p0("FILTER_GOPATH_ENRICH_MIN_EXPECTED", x)]])
        })
        observeEvent(input[[p0("FILTER_GOPATH_ENRICH_MIN_PRESENT", x)]], {
          req(input[[p0("FILTER_GOPATH_ENRICH_MIN_PRESENT", x)]])
          TabFilterGopathEnrichValues[[p0("FILTER_GOPATH_ENRICH_MIN_PRESENT_VALUE", x)]] <- 
            as.numeric(input[[p0("FILTER_GOPATH_ENRICH_MIN_PRESENT", x)]])
        })
        observeEvent(input[[p0("FILTER_GOPATH_ENRICH_MIN_TOTAL", x)]], {
          req(input[[p0("FILTER_GOPATH_ENRICH_MIN_TOTAL", x)]])
          TabFilterGopathEnrichValues[[p0("FILTER_GOPATH_ENRICH_MIN_TOTAL_VALUE", x)]] <- 
            as.numeric(input[[p0("FILTER_GOPATH_ENRICH_MIN_TOTAL", x)]])
        })
        observeEvent(input[[p0("FILTER_GOPATH_ENRICH_MAX_FP", x)]], {
          req(input[[p0("FILTER_GOPATH_ENRICH_MAX_FP", x)]])
          TabFilterGopathEnrichValues[[p0("FILTER_GOPATH_ENRICH_MAX_FP_VALUE", x)]] <- 
            as.numeric(input[[p0("FILTER_GOPATH_ENRICH_MAX_FP", x)]])
        })
        observeEvent(input[[p0("FILTER_GOPATH_ENRICH_MAX_FDR", x)]], {
          req(input[[p0("FILTER_GOPATH_ENRICH_MAX_FDR", x)]])
          TabFilterGopathEnrichValues[[p0("FILTER_GOPATH_ENRICH_MAX_FDR_VALUE", x)]] <- 
            as.numeric(input[[p0("FILTER_GOPATH_ENRICH_MAX_FDR", x)]])
        })
        if (!x %in% c("_KO", "_BLASTP_KO")) {
          observeEvent(input[[p0("FILTER_GOPATH_ENRICH_ONTOLOGY", x)]], {
            req(input[[p0("FILTER_GOPATH_ENRICH_ONTOLOGY", x)]])
            TabFilterGopathEnrichValues[[p0("FILTER_GOPATH_ENRICH_ONTOLOGY", x)]] <- 
              input[[p0("FILTER_GOPATH_ENRICH_ONTOLOGY", x)]]
          })
        }
      }
    )
  })

  #-----------------------------------------------------------------------------
  # Loading indicators

  TAB_HAS_LOADED <- reactiveValues(
    "KEGG Orthologs" = FALSE,
    "TIGRFAM" = FALSE,
    "PFAM" = FALSE,
    "Results - KO" = FALSE,
    "Results - TIGRFAM" = FALSE,
    "Results - PFAM" = FALSE,
    "BLASTP_INPUT_TAB" = TRUE
  )

  observeEvent(input$NAVBARPG, {
    req(input$NAVBARPG)
    if (!TAB_HAS_LOADED[[input$NAVBARPG]]) {
      # Sys.sleep(3)
      hide(anim = TRUE, animType = "fade", id = switch(input$NAVBARPG,
        "KEGG Orthologs" = "UI_KO_LOADING",
        "TIGRFAM" = "UI_TIGRFAM_LOADING",
        "PFAM" = "UI_PFAM_LOADING",
        "Results - KO" = "UI_BLASTP_KO_LOADING",
        "Results - TIGRFAM" = "UI_BLASTP_TIGRFAM_LOADING",
        "Results - PFAM" = "UI_BLASTP_PFAM_LOADING"
      ))
      TAB_HAS_LOADED[[input$NAVBARPG]] <- TRUE
    }
  })

}

msg("Finished loading server.R")
server
