# Independent pairwise community beta diversity. Base R + shiny only.
ptq_beta_metrics <- function(m, metric) {
  n <- nrow(m)
  out <- matrix(0,n,n,dimnames=list(rownames(m),rownames(m)))
  for (i in seq_len(n-1L)) for (j in seq.int(i+1L,n)) {
    x <- m[i,]; y <- m[j,]
    if (metric=='Bray-Curtis') {
      denom <- sum(x+y)
      val <- if(denom==0) NA_real_ else sum(abs(x-y))/denom
    } else {
      a <- sum(x>0 & y>0)
      b <- sum(x>0 & y==0)
      c <- sum(x==0 & y>0)
      if (metric=='Jaccard') val <- if(a+b+c==0) NA_real_ else (b+c)/(a+b+c)
      else if (metric=='Sørensen') val <- if(2*a+b+c==0) NA_real_ else (b+c)/(2*a+b+c)
      else if (metric=='Simpson turnover') val <- if(a+min(b,c)==0) NA_real_ else min(b,c)/(a+min(b,c))
      else stop('Unsupported beta diversity metric')
    }
    out[i,j] <- val; out[j,i] <- val
  }
  out
}
ptq_beta_ui <- function() {
  shiny::tabPanel('Community beta diversity',
    shiny::h3('Pairwise community beta diversity'),
    shiny::p('Upload a site-by-species abundance CSV. This module measures taxonomic or abundance-based beta diversity, not functional beta diversity. It requires no trait file.'),
    shiny::fluidRow(
      shiny::column(4,shiny::div(class='well',
        shiny::fileInput('beta_file','Upload site abundance CSV',accept='.csv'),
        shiny::selectInput('beta_metric','Dissimilarity index',
          choices=c('Jaccard','Sørensen','Simpson turnover','Bray-Curtis')),
        shiny::actionButton('beta_run','Calculate beta diversity',class='btn-primary'),
        shiny::downloadButton('beta_example','Download example abundance CSV'),
        shiny::p('Jaccard, Sørensen and Simpson turnover use presence/absence. Bray-Curtis uses nonnegative abundances.'),
        shiny::p('Pairs of two empty sites have undefined dissimilarity (NA). Ordination needs complete pairwise distances.'))),
      shiny::column(8,
        shiny::h4('Input checks'),shiny::verbatimTextOutput('beta_validation'),
        shiny::h4('Pairwise dissimilarity (first 12 sites)'),shiny::tableOutput('beta_matrix'),
        shiny::downloadButton('beta_matrix_csv','Download full dissimilarity matrix'),
        shiny::h4('Site PCoA ordination'),shiny::plotOutput('beta_pcoa',height='430px'),
        shiny::verbatimTextOutput('beta_notes'),
        shiny::h4('PCoA coordinates'),shiny::tableOutput('beta_scores'),
        shiny::downloadButton('beta_scores_csv','Download coordinates CSV'))))
}
ptq_beta_server <- function(input,output,session) {
  output$beta_example <- shiny::downloadHandler(
    filename=function()'site_abundance_example.csv',
    content=function(file)file.copy(file.path('data','site_abundance_example.csv'),file,overwrite=TRUE))
  data <- shiny::reactive({
    shiny::req(input$beta_file)
    d <- tryCatch(utils::read.csv(input$beta_file$datapath,check.names=FALSE,
      stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
      error=function(e)shiny::validate(shiny::need(FALSE,paste('Cannot read CSV:',conditionMessage(e)))))
    names(d) <- trimws(names(d))
    shiny::validate(shiny::need(!anyDuplicated(names(d)),'Duplicate column names.'),
      shiny::need('site_id' %in% names(d),'Required site_id column is missing.'),
      shiny::need(nrow(d)>=2,'At least two sites required.'),
      shiny::need(nrow(d)<=600,'Interactive pairwise analysis currently supports at most 600 sites.' ))
    d$site_id <- trimws(as.character(d$site_id))
    shiny::validate(shiny::need(all(!is.na(d$site_id) & nzchar(d$site_id)),'Missing site identifiers.'),
      shiny::need(!anyDuplicated(d$site_id),'Duplicate site identifiers.'))
    cols <- setdiff(names(d),'site_id')
    shiny::validate(shiny::need(length(cols)>0,'No species columns found.'))
    for(nm in cols) {
      v <- suppressWarnings(as.numeric(d[[nm]]))
      shiny::validate(shiny::need(all(is.na(d[[nm]]) | is.finite(v)),paste('Non-numeric abundance in:',nm)))
      d[[nm]] <- v
    }
    m <- as.matrix(d[cols]); storage.mode(m) <- 'double'; rownames(m) <- d$site_id
    shiny::validate(shiny::need(!anyNA(m),'Missing abundance cells: correct the CSV before analysis.'),
      shiny::need(all(is.finite(m) & m>=0),'Abundance must be finite and nonnegative.'))
    m
  })
  output$beta_validation <- shiny::renderText({
    if(is.null(input$beta_file))return('Upload an abundance CSV or download the example.')
    m <- data()
    paste0('Sites: ',nrow(m),' | Species columns: ',ncol(m),
      '\nEmpty sites: ',sum(rowSums(m)==0),
      '\nSpecies observed across dataset: ',sum(colSums(m)>0),
      '\nSelected metric: ',input$beta_metric,
      '\nAll abundances nonnegative; no missing cells.')
  })
  calc <- shiny::eventReactive(input$beta_run,{
    m <- data(); metric <- input$beta_metric
    dist <- ptq_beta_metrics(m,metric)
    list(m=m,dist=dist,metric=metric)
  })
  output$beta_matrix <- shiny::renderTable({
    d <- calc()$dist
    as.data.frame(round(d[seq_len(min(12,nrow(d))),seq_len(min(12,ncol(d))),drop=FALSE],4))
  },rownames=TRUE,striped=TRUE)
  output$beta_matrix_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_beta_dissimilarity.csv',
    content=function(file)utils::write.csv(calc()$dist,file,row.names=TRUE,na=''))
  ord <- shiny::reactive({
    z <- calc(); d <- z$dist
    shiny::validate(shiny::need(nrow(d)>=3,'At least three sites are required for two-axis PCoA.'),
      shiny::need(!anyNA(d),'Some site pairs have undefined dissimilarity (e.g. both empty). Fix the input before PCoA.'),
      shiny::need(any(d>0),'All site dissimilarities are zero. Ordination is undefined.'))
    fit <- tryCatch(stats::cmdscale(stats::as.dist(d),k=2,eig=TRUE,add=TRUE),
      error=function(e)shiny::validate(shiny::need(FALSE,paste('PCoA failed:',conditionMessage(e)))))
    shiny::validate(shiny::need(is.matrix(fit$points) && ncol(fit$points)>=2,
      'Two ordination axes could not be obtained.'))
    list(scores=data.frame(site_id=rownames(d),PCoA1=fit$points[,1],PCoA2=fit$points[,2]),fit=fit)
  })
  output$beta_pcoa <- shiny::renderPlot({
    x <- ord()$scores
    graphics::par(mar=c(4.5,4.5,2.5,1.5))
    graphics::plot(x$PCoA1,x$PCoA2,pch=19,col='#176d51',
      xlab='PCoA axis 1',ylab='PCoA axis 2',main='Community beta diversity ordination')
    graphics::abline(h=0,v=0,lty=3,col='grey70')
    if(nrow(x)<=40)graphics::text(x$PCoA1,x$PCoA2,labels=x$site_id,pos=3,cex=.7)
  },width=850,height=530,res=96)
  output$beta_notes <- shiny::renderText({
    z <- calc()
    paste0('Metric: ',z$metric,
      '\nDistance range: 0 to 1 for nonempty site comparisons.',
      '\nDissimilarity matrix is available even when PCoA cannot be calculated.',
      '\nPCoA uses Cailliez additive correction when needed.')
  })
  output$beta_scores <- shiny::renderTable(utils::head(ord()$scores,100),striped=TRUE)
  output$beta_scores_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_beta_PCoA_scores.csv',
    content=function(file)utils::write.csv(ord()$scores,file,row.names=FALSE))
}
