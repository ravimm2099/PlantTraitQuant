
# Community-weighted mean (CWM) trait centroid dissimilarity, an interpretable
# component of functional beta diversity; not a replacement for turnover partitioning.
ptq_fbeta_ui <- function() {
  shiny::tabPanel('Functional beta diversity',
    shiny::h3('Between-site functional differences in community trait centroids'),
    shiny::p('This first functional beta module calculates pairwise Euclidean distances between abundance-weighted community mean trait vectors (CWM). It is a centroid-based functional dissimilarity, not full functional turnover or nestedness.'),
    shiny::fluidRow(
      shiny::column(4,shiny::div(class='well',
        shiny::fileInput('fb_traits','Species traits CSV',accept='.csv'),
        shiny::fileInput('fb_abundance','Site abundance CSV',accept='.csv'),
        shiny::uiOutput('fb_trait_select'),
        shiny::checkboxInput('fb_presence','Use presence/absence weights instead of abundance',FALSE),
        shiny::actionButton('fb_run','Calculate functional differences',class='btn-primary'),
        shiny::downloadButton('fb_trait_example','Example traits'),
        shiny::downloadButton('fb_abundance_example','Example abundance'),
        shiny::p('Trait variables are standardized across species observed at any site. Sites with no species cannot be analysed.'))),
      shiny::column(8,
        shiny::h4('Validation'),shiny::verbatimTextOutput('fb_validation'),
        shiny::h4('Pairwise CWM distance matrix'),shiny::tableOutput('fb_distance_preview'),
        shiny::downloadButton('fb_distance_csv','Download distance matrix CSV'),
        shiny::h4('PCoA of community functional centroids'),
        shiny::plotOutput('fb_plot',height='430px'),
        shiny::h4('Site coordinates and community means'),
        shiny::tableOutput('fb_scores_preview'),
        shiny::downloadButton('fb_scores_csv','Download site results CSV'),
        shiny::h4('Method interpretation'),shiny::verbatimTextOutput('fb_notes'))))
}
ptq_fbeta_server <- function(input,output,session) {
  output$fb_trait_example <- shiny::downloadHandler(
    filename=function()'fd_traits_example.csv',
    content=function(file)file.copy(file.path('data','fd_traits_example.csv'),file,overwrite=TRUE))
  output$fb_abundance_example <- shiny::downloadHandler(
    filename=function()'fd_abundance_example.csv',
    content=function(file)file.copy(file.path('data','fd_abundance_example.csv'),file,overwrite=TRUE))
  readcsv <- function(upload) {
    shiny::req(upload)
    tryCatch(utils::read.csv(upload$datapath,check.names=FALSE,
      stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
      error=function(e)shiny::validate(shiny::need(FALSE,paste('CSV error:',conditionMessage(e)))))
  }
  datasets <- shiny::reactive({
    t <- readcsv(input$fb_traits)
    a <- readcsv(input$fb_abundance)
    names(t) <- trimws(names(t)); names(a) <- trimws(names(a))
    shiny::validate(
      shiny::need(!anyDuplicated(names(t)) && !anyDuplicated(names(a)),
        'Duplicate CSV column names.'),
      shiny::need('species_id' %in% names(t),'Traits require species_id.'),
      shiny::need('site_id' %in% names(a),'Abundance requires site_id.'),
      shiny::need(nrow(a)>=3,'At least three sites needed for two-axis ordination.'),
      shiny::need(nrow(a)<=500,'Limit of 500 sites for this interactive module.'))
    t$species_id <- trimws(as.character(t$species_id))
    a$site_id <- trimws(as.character(a$site_id))
    shiny::validate(
      shiny::need(all(!is.na(t$species_id) & nzchar(t$species_id)) &&
        !anyDuplicated(t$species_id),'Species IDs must be unique and nonempty.'),
      shiny::need(all(!is.na(a$site_id) & nzchar(a$site_id)) &&
        !anyDuplicated(a$site_id),'Site IDs must be unique and nonempty.'))
    ac <- setdiff(names(a),'site_id')
    shiny::validate(shiny::need(length(ac)>0,'Abundance file has no species columns.'))
    for(nm in ac) {
      raw <- a[[nm]]
      v <- suppressWarnings(as.numeric(raw))
      shiny::validate(shiny::need(all(is.na(raw) | is.finite(v)),
        paste('Invalid abundance values in',nm)))
      a[[nm]] <- v
    }
    m <- as.matrix(a[ac])
    storage.mode(m) <- 'double'
    shiny::validate(shiny::need(!anyNA(m),'Missing abundances must be corrected.'),
      shiny::need(all(is.finite(m) & m>=0),'Abundances must be finite and nonnegative.'),
      shiny::need(all(rowSums(m)>0),'Empty sites must be removed or populated.'))
    tc <- setdiff(names(t),'species_id')
    for(nm in tc) {
      raw <- t[[nm]]
      v <- suppressWarnings(as.numeric(raw))
      if(all(is.na(raw) | is.finite(v)))t[[nm]] <- v
    }
    numeric_cols <- tc[vapply(t[tc],is.numeric,logical(1))]
    shiny::validate(shiny::need(length(numeric_cols)>=2,'At least two numeric traits required.'))
    list(t=t,a=a,m=m,ac=ac,numeric_cols=numeric_cols)
  })
  output$fb_trait_select <- shiny::renderUI({
    x <- datasets()
    shiny::selectInput('fb_selected','Numeric traits to include',
      choices=x$numeric_cols,selected=x$numeric_cols,multiple=TRUE)
  })
  output$fb_validation <- shiny::renderText({
    if(is.null(input$fb_traits) || is.null(input$fb_abundance))
      return('Upload both files or download the example pair.')
    x <- datasets()
    used <- x$ac[colSums(x$m)>0]
    unmatched <- setdiff(used,x$t$species_id)
    chosen <- intersect(input$fb_selected,x$numeric_cols)
    paste0('Sites: ',nrow(x$a),'\nSpecies observed: ',length(used),
      '\nSelected numeric traits: ',length(chosen),
      '\nObserved species without traits: ',
      if(length(unmatched))paste(unmatched,collapse=', ') else 'none',
      '\nMissing selected trait values among all trait rows: ',
      if(length(chosen))sum(is.na(x$t[chosen])) else 0,
      '\nTrait records must be complete for every observed species.')
  })
  calculated <- shiny::eventReactive(input$fb_run,{
    x <- datasets()
    selected <- intersect(input$fb_selected,x$numeric_cols)
    shiny::validate(shiny::need(length(selected)>=2,'Select two or more numeric traits.'))
    used <- x$ac[colSums(x$m)>0]
    missing <- setdiff(used,x$t$species_id)
    shiny::validate(shiny::need(!length(missing),
      paste('Observed species absent from traits:',paste(missing,collapse=', '))))
    t <- x$t[match(used,x$t$species_id),selected,drop=FALSE]
    shiny::validate(shiny::need(all(stats::complete.cases(t)) &&
      all(is.finite(as.matrix(t))),'Selected traits have missing or nonfinite values for observed species.'))
    std <- vapply(t,stats::sd,numeric(1))
    shiny::validate(shiny::need(all(is.finite(std) & std>0),
      'Remove constant selected traits.'))
    traits <- scale(as.matrix(t))
    a <- x$m[,match(used,x$ac),drop=FALSE]
    if(isTRUE(input$fb_presence)) a <- (a>0)*1
    weights <- a/rowSums(a)
    cwm <- weights %*% traits
    rownames(cwm) <- x$a$site_id
    colnames(cwm) <- paste0('CWM_',selected)
    d <- as.matrix(stats::dist(cwm,method='euclidean'))
    ord <- stats::cmdscale(stats::as.dist(d),k=2,eig=TRUE)
    xy <- ord$points
    shiny::validate(shiny::need(!is.null(dim(xy)) && ncol(xy)>=2,
      'At least two nonzero ordination axes are required.'))
    result <- data.frame(site_id=x$a$site_id,PCoA1=xy[,1],PCoA2=xy[,2],
      cwm,check.names=FALSE)
    list(d=d,result=result,sites=nrow(a),species=length(used),
      traits=length(selected))
  })
  output$fb_distance_preview <- shiny::renderTable({
    d <- calculated()$d
    as.data.frame(d[seq_len(min(10,nrow(d))),seq_len(min(10,ncol(d))),drop=FALSE])
  },rownames=TRUE,striped=TRUE)
  output$fb_scores_preview <- shiny::renderTable(utils::head(calculated()$result,100),striped=TRUE)
  output$fb_plot <- shiny::renderPlot({
    x <- calculated()$result
    graphics::par(mar=c(4.5,4.5,2.5,1.5))
    graphics::plot(x$PCoA1,x$PCoA2,pch=19,col='#16735b',cex=1.2,
      xlab='PCoA axis 1',ylab='PCoA axis 2',
      main='Differences between community trait centroids')
    graphics::abline(h=0,v=0,lty=3,col='grey70')
    if(nrow(x)<=50)graphics::text(x$PCoA1,x$PCoA2,labels=x$site_id,pos=3,cex=.7)
  },width=850,height=520,res=96)
  output$fb_distance_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_CWM_functional_distances.csv',
    content=function(file)utils::write.csv(calculated()$d,file,row.names=TRUE))
  output$fb_scores_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_CWM_site_scores.csv',
    content=function(file)utils::write.csv(calculated()$result,file,row.names=FALSE))
  output$fb_notes <- shiny::renderText({
    z <- calculated()
    paste0('Method: Euclidean distance between community-weighted means of ',
      z$traits,' globally standardized numeric traits.\n',
      'Sites: ',z$sites,' | Observed species: ',z$species,
      '\nWeighting: ',if(isTRUE(input$fb_presence))'presence/absence' else 'relative abundance',
      '\nThis captures centroid shifts but not within-community dispersion,',
      ' functional turnover partitioning, or nestedness.\n',
      'Zero distances indicate equal centroids, not necessarily identical species composition.')
  })
}
