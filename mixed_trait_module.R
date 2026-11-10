
# Independent mixed-trait Gower dissimilarity and PCoA
# Only base R and shiny; no external ecology packages.
ptq_gower_distance <- function(d, trait_cols) {
  n <- nrow(d)
  sums <- matrix(0, n, n)
  weights <- matrix(0, n, n)
  dropped <- character()
  for (nm in trait_cols) {
    v <- d[[nm]]
    if (is.numeric(v)) {
      v[!is.finite(v)] <- NA_real_
      observed <- v[!is.na(v)]
      span <- if(length(observed)) diff(range(observed)) else 0
      if (!is.finite(span) || span <= 0) { dropped <- c(dropped,nm); next }
      delta <- abs(outer(v,v,"-")) / span
    } else {
      v <- trimws(as.character(v))
      v[is.na(v) | !nzchar(v)] <- NA_character_
      if (length(unique(v[!is.na(v)])) < 2L) { dropped <- c(dropped,nm); next }
      delta <- outer(v,v,"!=") * 1
    }
    ok <- outer(!is.na(v),!is.na(v),"&")
    delta[!ok] <- 0
    sums <- sums + delta
    weights <- weights + ok
  }
  if (!any(weights[upper.tri(weights)]>0))
    stop("No comparable trait pairs. Check missing values and selected traits.")
  dist_mat <- sums / pmax(weights,1)
  dist_mat[weights==0] <- NA_real_
  diag(dist_mat) <- 0
  dimnames(dist_mat) <- list(d$species_id,d$species_id)
  list(distance=dist_mat,excluded_traits=dropped,
       uncomparable=sum(is.na(dist_mat[upper.tri(dist_mat)])))
}
ptq_mixed_ui <- function() {
  shiny::tabPanel('Mixed traits (Gower + PCoA)',
    shiny::h3('Mixed numeric and categorical plant traits'),
    shiny::p('Independent analysis: upload one species-trait CSV, choose traits and calculate Gower distances and PCoA. No abundance file required.'),
    shiny::fluidRow(
      shiny::column(4,shiny::div(class='well',
        shiny::fileInput('mixed_file','Upload mixed species traits CSV',accept='.csv'),
        shiny::uiOutput('mixed_trait_ui'),
        shiny::actionButton('mixed_run','Calculate Gower + PCoA',class='btn-primary'),
        shiny::downloadButton('mixed_example','Download mixed-trait example'),
        shiny::p('Numeric traits are scaled by their observed range. Categorical traits use 0/1 mismatch. Missing values are handled pairwise.'),
        shiny::p('Ordination needs at least three species and complete pairwise distances.'))),
      shiny::column(8,
        shiny::h4('Input and missing-value checks'),shiny::verbatimTextOutput('mixed_checks'),
        shiny::h4('PCoA ordination'),shiny::plotOutput('mixed_plot',height='450px'),
        shiny::h4('Ordination diagnostics'),shiny::verbatimTextOutput('mixed_diagnostics'),
        shiny::h4('Species coordinates'),shiny::tableOutput('mixed_scores'),
        shiny::downloadButton('mixed_scores_csv','Download coordinates CSV'),
        shiny::h4('Gower distance matrix (preview)'),shiny::tableOutput('mixed_distances'),
        shiny::downloadButton('mixed_distance_csv','Download distance matrix CSV'))))
}
ptq_mixed_server <- function(input,output,session) {
  example <- file.path('data','mixed_traits_example.csv')
  output$mixed_example <- shiny::downloadHandler(
    filename=function() 'mixed_traits_example.csv',
    content=function(file) file.copy(example,file,overwrite=TRUE))
  dat <- shiny::reactive({
    shiny::req(input$mixed_file)
    d <- tryCatch(utils::read.csv(input$mixed_file$datapath,check.names=FALSE,
                    stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),
                    fileEncoding='UTF-8-BOM'),
      error=function(e) shiny::validate(shiny::need(FALSE,paste('Cannot read CSV:',conditionMessage(e)))))
    names(d) <- trimws(names(d))
    shiny::validate(shiny::need('species_id' %in% names(d),'A species_id column is required.'),
      shiny::need(!anyDuplicated(names(d)),'Duplicate column names found.'),
      shiny::need(nrow(d)>=3,'At least three species are required.'))
    d$species_id <- trimws(as.character(d$species_id))
    shiny::validate(shiny::need(all(!is.na(d$species_id) & nzchar(d$species_id)),
      'Missing species identifiers.'),shiny::need(!anyDuplicated(d$species_id),
      'Duplicate species identifiers.'))
    cols <- setdiff(names(d),'species_id')
    shiny::validate(shiny::need(length(cols)>0,'Add at least one trait column.'))
    for(nm in cols) {
      if(is.character(d[[nm]])) {
        v <- suppressWarnings(as.numeric(d[[nm]]))
        nonempty <- !is.na(d[[nm]]) & nzchar(trimws(d[[nm]]))
        if(all(!nonempty | is.finite(v))) d[[nm]] <- v
      }
    }
    d
  })
  output$mixed_trait_ui <- shiny::renderUI({
    d <- dat()
    cols <- setdiff(names(d),'species_id')
    shiny::selectInput('mixed_traits','Choose traits',choices=cols,selected=cols,multiple=TRUE)
  })
  output$mixed_checks <- shiny::renderText({
    if(is.null(input$mixed_file)) return('Upload a mixed-trait CSV or download the example.')
    d <- dat()
    cols <- intersect(input$mixed_traits,setdiff(names(d),'species_id'))
    paste0('Species: ',nrow(d),'\nSelected traits: ',length(cols),
      '\nNumeric: ',sum(vapply(d[cols],is.numeric,logical(1))),
      ' | Categorical: ',sum(!vapply(d[cols],is.numeric,logical(1))),
      '\nMissing selected cells: ',sum(is.na(d[cols])),
      '\nRows with at least one missing selected trait: ',sum(!stats::complete.cases(d[cols])),
      '\nNo fixed number of trait columns; pairwise distances require comparable values.')
  })
  result <- shiny::eventReactive(input$mixed_run,{
    d <- dat()
    cols <- intersect(input$mixed_traits,setdiff(names(d),'species_id'))
    shiny::validate(shiny::need(length(cols)>0,'Select one or more traits.'),
      shiny::need(nrow(d)<=1200,'For interactive Gower + PCoA, use at most 1200 species per run. This is a computational safeguard, not a limit on other tabs.'))
    z <- tryCatch(ptq_gower_distance(d,cols),
      error=function(e) shiny::validate(shiny::need(FALSE,conditionMessage(e))))
    shiny::validate(shiny::need(z$uncomparable==0,
      paste('There are',z$uncomparable,'species pairs with no shared measured traits. Supply at least one comparable trait per pair before PCoA.')))
    ord <- tryCatch(stats::cmdscale(stats::as.dist(z$distance),k=2,eig=TRUE,add=TRUE),
      error=function(e) shiny::validate(shiny::need(FALSE,paste('PCoA failed:',conditionMessage(e)))))
    coords <- ord$points
    shiny::validate(shiny::need(!is.null(dim(coords)) && ncol(coords)>=2,
      'PCoA did not produce two axes. Add informative traits/species.'))
    pos <- ord$eig[ord$eig>1e-10]
    variance <- if(length(pos)) 100*ord$eig[1:2]/sum(pos) else c(NA_real_,NA_real_)
    scores <- data.frame(species_id=d$species_id,PCoA1=coords[,1],PCoA2=coords[,2])
    list(gower=z,ord=ord,scores=scores,variance=variance)
  })
  output$mixed_diagnostics <- shiny::renderText({
    z <- result()
    paste0('Species: ',nrow(z$scores),
      '\nTraits excluded as constant/all-missing: ',
      if(length(z$gower$excluded_traits)) paste(z$gower$excluded_traits,collapse=', ') else 'none',
      '\nPairs lacking shared observed traits: ',z$gower$uncomparable,
      '\nPCoA uses Cailliez correction when required (stats::cmdscale(add=TRUE)).',
      '\nPCoA axis 1: ',round(z$variance[1],2),'% of positive eigenvalue sum',
      '\nPCoA axis 2: ',round(z$variance[2],2),'% of positive eigenvalue sum')
  })
  output$mixed_plot <- shiny::renderPlot({
    z <- result()
    x <- z$scores
    graphics::par(mar=c(4.5,4.5,2,1.5))
    graphics::plot(x$PCoA1,x$PCoA2,pch=19,col='#146c51',cex=1.1,
      xlab='PCoA axis 1',ylab='PCoA axis 2',main='Gower trait-space ordination')
    graphics::abline(h=0,v=0,lty=3,col='grey70')
    if(nrow(x)<=50) graphics::text(x$PCoA1,x$PCoA2,labels=x$species_id,pos=3,cex=.65)
  },width=850,height=530,res=96)
  output$mixed_scores <- shiny::renderTable(utils::head(result()$scores,100),striped=TRUE)
  output$mixed_distances <- shiny::renderTable({
    m <- result()$gower$distance
    as.data.frame(m[seq_len(min(8,nrow(m))),seq_len(min(8,ncol(m))),drop=FALSE])
  },rownames=TRUE,striped=TRUE)
  output$mixed_scores_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_Gower_PCoA_scores.csv',
    content=function(file)utils::write.csv(result()$scores,file,row.names=FALSE))
  output$mixed_distance_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_Gower_distance.csv',
    content=function(file)utils::write.csv(result()$gower$distance,file,row.names=TRUE))
}
