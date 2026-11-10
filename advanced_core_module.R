
# Advanced germplasm collection: PAM-style medoid swap optimization and strata quotas.
# Depends on ptq_core_gower from core_collection_module.R; no extra R packages.
ptq_medoid_cost <- function(d, ids) mean(apply(d[,ids,drop=FALSE],1,min))
ptq_medoid_opt <- function(d,k,starts=3L,seed=123L,max_pass=10L) {
  n <- nrow(d)
  if(k>=n) return(seq_len(n))
  best <- NULL; best_cost <- Inf
  for(r in seq_len(starts)) {
    set.seed(seed+r-1L)
    med <- sample.int(n,k)
    cost <- ptq_medoid_cost(d,med)
    for(pass in seq_len(max_pass)) {
      improved <- FALSE
      # PAM-like single-medoid swap local search, not cluster::pam().
      for(pos in seq_len(k)) {
        candidates <- setdiff(seq_len(n),med)
        for(cand in candidates) {
          trial <- med; trial[pos] <- cand
          score <- ptq_medoid_cost(d,trial)
          if(score < cost-1e-10) {
            med <- trial;cost <- score;improved <- TRUE
          }
        }
      }
      if(!improved)break
    }
    if(cost<best_cost) {best <- med;best_cost <- cost}
  }
  best
}
ptq_allocate_strata <- function(group,k) {
  tab <- table(group); n <- sum(tab)
  quota <- floor(as.numeric(tab)*k/n)
  quota <- pmin(quota,as.integer(tab))
  remaining <- k-sum(quota)
  remainder <- as.numeric(tab)*k/n-quota
  while(remaining>0) {
    candidates <- which(quota<as.integer(tab))
    j <- candidates[which.max(remainder[candidates])]
    quota[j] <- quota[j]+1L
    remainder[j] <- -1
    remaining <- remaining-1L
  }
  setNames(as.integer(quota),names(tab))
}
ptq_advanced_core_ui <- function() {
  shiny::tabPanel('Advanced core selection',
    shiny::h3('Core and mini-core: medoid optimization and stratification'),
    shiny::p('Compare representative germplasm subsets using Gower descriptors, multiple-start medoid swaps and optional proportional strata.'),
    shiny::fluidRow(
      shiny::column(4,shiny::div(class='well',
        shiny::fileInput('ac_file','Accession descriptor CSV',accept='.csv'),
        shiny::uiOutput('ac_group_ui'),
        shiny::uiOutput('ac_traits_ui'),
        shiny::numericInput('ac_core','Core size (%)',20,min=1,max=100),
        shiny::numericInput('ac_mini','Mini-core size (%)',5,min=1,max=100),
        shiny::numericInput('ac_starts','Random starts',3,min=1,max=10,step=1),
        shiny::numericInput('ac_seed','Random seed',123,min=1,max=1000000),
        shiny::actionButton('ac_go','Optimize selections',class='btn-primary'),
        shiny::downloadButton('ac_example','Download example accessions'),
        shiny::p('Mini-core is optimized within the core. Stratified core selection allocates quotas proportional to group size.'))),
      shiny::column(8,
        shiny::h4('Data and design checks'),shiny::verbatimTextOutput('ac_checks'),
        shiny::h4('Collection coverage and retention'),shiny::tableOutput('ac_summary'),
        shiny::h4('Accession selection'),shiny::tableOutput('ac_selected'),
        shiny::downloadButton('ac_export','Download selection CSV'),
        shiny::h4('Descriptor retention by collection'),shiny::tableOutput('ac_retention'),
        shiny::downloadButton('ac_retention_csv','Download descriptor retention CSV'),
        shiny::h4('Method limitations'),shiny::verbatimTextOutput('ac_notes'))))
}
ptq_advanced_core_server <- function(input,output,session) {
  output$ac_example <- shiny::downloadHandler(
    filename=function()'germplasm_accessions_example.csv',
    content=function(file)file.copy(file.path('data','germplasm_accessions_example.csv'),file,overwrite=TRUE))
  data <- shiny::reactive({
    shiny::req(input$ac_file)
    d <- tryCatch(utils::read.csv(input$ac_file$datapath,check.names=FALSE,
      stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
      error=function(e)shiny::validate(shiny::need(FALSE,paste('CSV error:',conditionMessage(e)))))
    names(d) <- trimws(names(d))
    shiny::validate(shiny::need('accession_id' %in% names(d),'Missing accession_id.'),
      shiny::need(!anyDuplicated(names(d)),'Duplicate column names.'),
      shiny::need(nrow(d)>=5 && nrow(d)<=350,
        'Advanced optimization accepts 5 to 350 accessions. Larger collections need a faster solver.'))
    d$accession_id <- trimws(as.character(d$accession_id))
    shiny::validate(shiny::need(!anyNA(d$accession_id) && all(nzchar(d$accession_id)) &&
      !anyDuplicated(d$accession_id),'Accession IDs must be unique and nonempty.'))
    cols <- setdiff(names(d),'accession_id')
    for(nm in cols)if(is.character(d[[nm]])) {
      raw <- d[[nm]]; v <- suppressWarnings(as.numeric(raw))
      nonempty <- !is.na(raw) & nzchar(trimws(raw))
      if(all(!nonempty | is.finite(v))) d[[nm]] <- v
    }
    list(d=d,cols=cols)
  })
  output$ac_group_ui <- shiny::renderUI({
    x <- data()
    shiny::selectInput('ac_group','Optional stratum/group column',
      choices=c('No stratification'='',x$cols),selected='')
  })
  output$ac_traits_ui <- shiny::renderUI({
    x <- data()
    cols <- setdiff(x$cols,input$ac_group)
    shiny::selectInput('ac_traits','Gower descriptors',choices=cols,
      selected=cols,multiple=TRUE)
  })
  output$ac_checks <- shiny::renderText({
    if(is.null(input$ac_file))return('Upload an accession CSV or download the example.')
    x <- data()
    cols <- intersect(input$ac_traits,setdiff(x$cols,input$ac_group))
    paste0('Accessions: ',nrow(x$d),'\nSelected descriptors: ',length(cols),
      '\nMissing descriptor cells: ',sum(is.na(x$d[cols])),
      '\nStratification: ',if(nzchar(input$ac_group))input$ac_group else 'none',
      '\nSelection uses equal-weight Gower distances; all accession pairs must be comparable.')
  })
  calc <- shiny::eventReactive(input$ac_go,{
    x <- data();d <- x$d;n <- nrow(d)
    cols <- intersect(input$ac_traits,setdiff(x$cols,input$ac_group))
    shiny::validate(shiny::need(length(cols)>0,'Select descriptors.'),
      shiny::need(is.finite(input$ac_core) && input$ac_core>0 && input$ac_core<=100,'Invalid core percentage.'),
      shiny::need(is.finite(input$ac_mini) && input$ac_mini>0 &&
        input$ac_mini<=input$ac_core,'Mini-core percentage must be positive and no larger than core percentage.'))
    kc <- max(1L,min(n,ceiling(n*input$ac_core/100)))
    km <- max(1L,min(kc,ceiling(n*input$ac_mini/100)))
    starts <- as.integer(input$ac_starts)
    seed <- as.integer(input$ac_seed)
    shiny::validate(shiny::need(!is.na(starts) && starts>=1 && starts<=10,'Starts must be 1 to 10.'),
      shiny::need(!is.na(seed) && seed>=1,'Seed must be a positive integer.'))
    g <- tryCatch(ptq_core_gower(d,cols),
      error=function(e)shiny::validate(shiny::need(FALSE,conditionMessage(e))))
    shiny::validate(shiny::need(g$missing_pairs==0,
      paste(g$missing_pairs,'accession pairs have no comparable traits.')))
    dist <- g$distance
    if(!is.null(input$ac_group) && nzchar(input$ac_group)) {
      group <- trimws(as.character(d[[input$ac_group]]))
      shiny::validate(shiny::need(all(!is.na(group) & nzchar(group)),
        'Group labels cannot be missing or blank.'))
      quota <- ptq_allocate_strata(group,kc)
      core <- integer()
      for(label in names(quota)) {
        idx <- which(group==label)
        if(quota[[label]]>0) {
          selected <- ptq_medoid_opt(dist[idx,idx,drop=FALSE],
            quota[[label]],starts,seed)
          core <- c(core,idx[selected])
        }
      }
    } else {
      core <- ptq_medoid_opt(dist,kc,starts,seed)
    }
    mini_local <- ptq_medoid_opt(dist[core,core,drop=FALSE],km,starts,seed+100L)
    mini <- core[mini_local]
    nc <- apply(dist[,core,drop=FALSE],1,min)
    nm <- apply(dist[,mini,drop=FALSE],1,min)
    selection <- data.frame(accession_id=d$accession_id,
      collection=ifelse(seq_len(n)%in%mini,'Mini-core (also core)',
        ifelse(seq_len(n)%in%core,'Core only','Not selected')),
      stringsAsFactors=FALSE)
    summaries <- data.frame(collection=c('Core','Mini-core'),
      selected=c(kc,km),
      percent=c(100*kc/n,100*km/n),
      mean_nearest_Gower=c(mean(nc),mean(nm)),
      max_nearest_Gower=c(max(nc),max(nm)))
    retention <- list()
    for(label in c('Core','Mini-core')) {
      idx <- if(label=='Core')core else mini
      for(nm0 in cols) {
        v <- d[[nm0]]
        if(is.numeric(v)) {
          original <- range(v,na.rm=TRUE)
          subset <- range(v[idx],na.rm=TRUE)
          score <- if(all(is.finite(c(original,subset))) &&
             diff(original)>0) 100*diff(subset)/diff(original) else NA_real_
          measure <- 'Observed numeric range retention (%)'
        } else {
          orig <- unique(as.character(v[!is.na(v)]))
          sub <- unique(as.character(v[idx][!is.na(v[idx])]))
          score <- if(length(orig))100*length(intersect(orig,sub))/length(orig) else NA_real_
          measure <- 'Observed category retention (%)'
        }
        retention[[length(retention)+1L]] <- data.frame(
          collection=label,descriptor=nm0,measure=measure,retention_percent=score)
      }
    }
    list(selection=selection,summary=summaries,
      retention=do.call(rbind,retention),excluded=g$omitted,
      group=if(nzchar(input$ac_group))input$ac_group else 'none',
      starts=starts)
  })
  output$ac_summary <- shiny::renderTable(calc()$summary,striped=TRUE,digits=4)
  output$ac_selected <- shiny::renderTable(utils::head(
    calc()$selection[calc()$selection$collection!='Not selected',,drop=FALSE],100),
    striped=TRUE)
  output$ac_retention <- shiny::renderTable(calc()$retention,striped=TRUE,digits=3)
  output$ac_export <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_advanced_core_selection.csv',
    content=function(file)utils::write.csv(calc()$selection,file,row.names=FALSE))
  output$ac_retention_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_descriptor_retention.csv',
    content=function(file)utils::write.csv(calc()$retention,file,row.names=FALSE,na=''))
  output$ac_notes <- shiny::renderText({
    z <- calc()
    paste0('Algorithm: multiple-start single-medoid-swap local search (PAM-style, not cluster::pam).\n',
      'Objective: minimize mean nearest-representative Gower distance.\n',
      'Random starts: ',z$starts,' | Stratification: ',z$group,
      '\nMini-core is optimized from the core and is always its subset.\n',
      'Retention measures observed descriptor ranges or categorical levels, not allele retention.\n',
      'A group quota of zero is possible for small strata; use dedicated quotas if every group must be represented.\n',
      'No global optimum guarantee or independent reference benchmark. Excluded traits: ',
      if(length(z$excluded))paste(z$excluded,collapse=', ') else 'none')
  })
}
