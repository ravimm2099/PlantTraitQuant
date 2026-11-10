# Core and mini-core collection using deterministic greedy k-center (maximin)
# No additional packages. This is a diversity-coverage heuristic, not PAM.
ptq_core_ui <- function() {
  shiny::tabPanel('Core and mini-core collection',
    shiny::h3('Germplasm core and mini-core selection'),
    shiny::p('Select representative accessions from numerical and categorical descriptors using Gower dissimilarity and greedy maximum-distance coverage.'),
    shiny::fluidRow(
      shiny::column(4,shiny::div(class='well',
        shiny::fileInput('core_file','Upload accession descriptor CSV',accept='.csv'),
        shiny::uiOutput('core_trait_ui'),
        shiny::numericInput('core_pct','Core size (% of accessions)',20,min=1,max=100,step=1),
        shiny::numericInput('mini_pct','Mini-core size (% of accessions)',5,min=1,max=100,step=1),
        shiny::numericInput('core_seed','Random seed for first accession',123,min=1,max=1000000,step=1),
        shiny::actionButton('core_go','Select core and mini-core',class='btn-primary'),
        shiny::downloadButton('core_example','Download example accessions'),
        shiny::p('Mini-core is a subset of the core. No population grouping or geographic stratification is assumed.'))),
      shiny::column(8,
        shiny::h4('Data checks'),shiny::verbatimTextOutput('core_checks'),
        shiny::h4('Selection summary'),shiny::tableOutput('core_summary'),
        shiny::h4('Selected accessions'),shiny::tableOutput('core_preview'),
        shiny::downloadButton('core_export','Download selected accession IDs'),
        shiny::downloadButton('core_full_export','Download selected descriptor records'),
        shiny::h4('Coverage diagnostics'),shiny::plotOutput('core_coverage_plot',height='380px'),
        shiny::downloadButton('core_dist_export','Download nearest-representative distances'),
        shiny::h4('Method notes'),shiny::verbatimTextOutput('core_notes'))))
}
ptq_core_gower <- function(d,cols) {
  n <- nrow(d)
  sums <- matrix(0,n,n); weights <- matrix(0,n,n)
  omitted <- character()
  for(nm in cols) {
    v <- d[[nm]]
    if(is.numeric(v)) {
      v[!is.finite(v)] <- NA_real_
      good <- v[!is.na(v)]
      span <- if(length(good))diff(range(good)) else 0
      if(!is.finite(span) || span<=0) {omitted <- c(omitted,nm);next}
      delta <- abs(outer(v,v,'-'))/span
    } else {
      v <- trimws(as.character(v));v[is.na(v)|!nzchar(v)] <- NA_character_
      if(length(unique(v[!is.na(v)]))<2) {omitted <- c(omitted,nm);next}
      delta <- outer(v,v,'!=')*1
    }
    valid <- outer(!is.na(v),!is.na(v),'&')
    delta[!valid] <- 0
    sums <- sums+delta;weights <- weights+valid
  }
  if(!any(weights[upper.tri(weights)]>0))stop('No informative comparable traits.')
  out <- sums/pmax(weights,1);out[weights==0] <- NA_real_
  diag(out) <- 0
  list(distance=out,omitted=omitted,missing_pairs=sum(is.na(out[upper.tri(out)])))
}
ptq_core_maximin <- function(distmat,k,initial=NULL,seed=123L) {
  n <- nrow(distmat)
  if(is.null(initial)) {
    set.seed(seed);initial <- sample.int(n,1)
  }
  chosen <- as.integer(initial)
  while(length(chosen)<k) {
    nearest <- apply(distmat[,chosen,drop=FALSE],1,min)
    nearest[chosen] <- -Inf
    chosen <- c(chosen,which.max(nearest))
  }
  chosen
}
ptq_core_server <- function(input,output,session) {
  output$core_example <- shiny::downloadHandler(
    filename=function()'germplasm_accessions_example.csv',
    content=function(file)file.copy(file.path('data','germplasm_accessions_example.csv'),file,overwrite=TRUE))
  data <- shiny::reactive({
    shiny::req(input$core_file)
    d <- tryCatch(utils::read.csv(input$core_file$datapath,check.names=FALSE,
       stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
       error=function(e)shiny::validate(shiny::need(FALSE,paste('CSV error:',conditionMessage(e)))))
    names(d) <- trimws(names(d))
    shiny::validate(shiny::need('accession_id' %in% names(d),'Required column accession_id is missing.'),
      shiny::need(!anyDuplicated(names(d)),'Duplicate column names.'),
      shiny::need(nrow(d)>=5,'At least five accessions required.'),
      shiny::need(nrow(d)<=1500,'Interactive limit is 1500 accessions (quadratic distance memory).'))
    d$accession_id <- trimws(as.character(d$accession_id))
    shiny::validate(shiny::need(all(!is.na(d$accession_id)&nzchar(d$accession_id)) &&
      !anyDuplicated(d$accession_id),'Accession IDs must be unique and nonempty.'))
    cols <- setdiff(names(d),'accession_id')
    shiny::validate(shiny::need(length(cols)>0,'Add at least one descriptor column.'))
    for(nm in cols) if(is.character(d[[nm]])) {
      raw <- d[[nm]];v <- suppressWarnings(as.numeric(raw))
      nonempty <- !is.na(raw)&nzchar(trimws(raw))
      if(all(!nonempty|is.finite(v)))d[[nm]] <- v
    }
    list(d=d,cols=cols)
  })
  output$core_trait_ui <- shiny::renderUI({
    x <- data()
    shiny::selectInput('core_cols','Select accession descriptors',choices=x$cols,
      selected=x$cols,multiple=TRUE)
  })
  output$core_checks <- shiny::renderText({
    if(is.null(input$core_file))return('Upload an accession CSV or use the example.')
    x <- data();cols <- intersect(input$core_cols,x$cols)
    paste0('Accessions: ',nrow(x$d),'\nSelected descriptors: ',length(cols),
      '\nNumeric descriptors: ',sum(vapply(x$d[cols],is.numeric,logical(1))),
      '\nCategorical descriptors: ',sum(!vapply(x$d[cols],is.numeric,logical(1))),
      '\nMissing selected values: ',sum(is.na(x$d[cols])),
      '\nPairwise missing descriptors are ignored, but every pair must share an observed informative descriptor.')
  })
  calc <- shiny::eventReactive(input$core_go,{
    x <- data();cols <- intersect(input$core_cols,x$cols)
    shiny::validate(shiny::need(length(cols)>0,'Select at least one descriptor.'),
      shiny::need(is.finite(input$core_pct) && input$core_pct>0 && input$core_pct<=100,'Core percentage must be 1 to 100.'),
      shiny::need(is.finite(input$mini_pct) && input$mini_pct>0 && input$mini_pct<=100,'Mini-core percentage must be 1 to 100.'),
      shiny::need(input$mini_pct<=input$core_pct,'Mini-core percentage must not exceed core percentage.'))
    g <- tryCatch(ptq_core_gower(x$d,cols),error=function(e)shiny::validate(shiny::need(FALSE,conditionMessage(e))))
    shiny::validate(shiny::need(g$missing_pairs==0,
      paste(g$missing_pairs,'accession pairs share no observed informative descriptors. Correct missing values before selection.')))
    n <- nrow(x$d);kc <- min(n,max(1,ceiling(n*input$core_pct/100)))
    km <- min(kc,max(1,ceiling(n*input$mini_pct/100)))
    seed <- as.integer(input$core_seed)
    shiny::validate(shiny::need(!is.na(seed),'Seed must be a valid integer.'))
    order <- ptq_core_maximin(g$distance,kc,seed=seed)
    # Nested mini-core is the prefix of the deterministic core-selection sequence.
    mini <- order[seq_len(km)]
    near_core <- apply(g$distance[,order,drop=FALSE],1,min)
    near_mini <- apply(g$distance[,mini,drop=FALSE],1,min)
    group <- rep('Not selected',n)
    group[order] <- 'Core only';group[mini] <- 'Mini-core (also core)'
    selection <- data.frame(accession_id=x$d$accession_id,
      selection=group,core_rank=match(seq_len(n),order),
      mini_core_rank=match(seq_len(n),mini),stringsAsFactors=FALSE)
    coverage <- data.frame(accession_id=x$d$accession_id,
      nearest_core_distance=near_core,nearest_mini_core_distance=near_mini)
    summary <- data.frame(collection=c('Core','Mini-core'),
      selected=c(kc,km),percent_of_total=round(100*c(kc,km)/n,2),
      mean_nearest_distance=c(mean(near_core),mean(near_mini)),
      max_nearest_distance=c(max(near_core),max(near_mini)))
    list(selection=selection,coverage=coverage,summary=summary,
      selected_records=cbind(selection[selection$selection!='Not selected',,drop=FALSE],
        x$d[selection$selection!='Not selected',setdiff(names(x$d),'accession_id'),drop=FALSE]),
      excluded=g$omitted,n=n,core=kc,mini=km)
  })
  output$core_summary <- shiny::renderTable(calc()$summary,striped=TRUE,digits=4)
  output$core_preview <- shiny::renderTable(utils::head(calc()$selection[calc()$selection!='Not selected',],100),striped=TRUE)
  output$core_coverage_plot <- shiny::renderPlot({
    z <- calc();a <- z$coverage
    graphics::par(mar=c(4.5,4.5,2.5,1.5))
    graphics::boxplot(a$nearest_core_distance,a$nearest_mini_core_distance,
      names=c('Core','Mini-core'),ylab='Nearest representative Gower distance',
      main='Descriptor coverage of original accessions',col=c('#a4d9bf','#e9cba1'))
  },width=850,height=450,res=96)
  output$core_export <- shiny::downloadHandler(filename=function()'PlantTraitQuant_core_selection.csv',
    content=function(file)utils::write.csv(calc()$selection,file,row.names=FALSE,na=''))
  output$core_full_export <- shiny::downloadHandler(filename=function()'PlantTraitQuant_selected_accessions.csv',
    content=function(file)utils::write.csv(calc()$selected_records,file,row.names=FALSE,na=''))
  output$core_dist_export <- shiny::downloadHandler(filename=function()'PlantTraitQuant_core_coverage.csv',
    content=function(file)utils::write.csv(calc()$coverage,file,row.names=FALSE))
  output$core_notes <- shiny::renderText({
    z <- calc()
    paste0('Selection: seeded greedy maximum-distance (k-center) heuristic. Not PAM or optimization of average dissimilarity.\n',
      'Core: ',z$core,' of ',z$n,' accessions; mini-core: ',z$mini,' (subset of core).\n',
      'Each descriptor contributes equal weight to pairwise Gower distance. Numeric differences are divided by observed range; categorical differences are 0/1.\n',
      'Excluded constant/all-missing descriptors: ',if(length(z$excluded))paste(z$excluded,collapse=', ') else 'none',
      '\nCoverage summaries measure nearest-selected descriptor distance, not allele retention or geographic representation.\n',
      'This version does not implement stratified sampling, PAM, genetic-marker diversity, or external benchmark validation.')
  })
}
