# Two-dimensional FDiv following the convex-hull-vertex-centroid approach.
ptq_fdiv_value <- function(x, abundance) {
  ids <- which(abundance > 0)
  if(length(ids) < 3L) return(NA_real_)
  xy <- x[ids,,drop=FALSE]
  if(qr(cbind(1,xy))$rank < 3L) return(NA_real_)
  hull <- grDevices::chull(xy[,1],xy[,2])
  centre <- colMeans(xy[hull,,drop=FALSE])
  radii <- sqrt(rowSums(sweep(xy,2,centre,'-')^2))
  dbar <- mean(radii)
  p <- abundance[ids]/sum(abundance[ids])
  delta <- sum(p*(radii-dbar))
  delta_abs <- sum(p*abs(radii-dbar))
  den <- dbar+delta_abs
  if(!is.finite(den) || den<=0) return(NA_real_)
  max(0,min(1,(dbar+delta)/den))
}
ptq_fdiv_ui <- function() {
  shiny::tabPanel('Functional divergence (2D FDiv)',
    shiny::h3('Functional divergence (FDiv) in 2D PCA trait space'),
    shiny::p('Independent linked trait and abundance files. Calculates an abundance-weighted divergence relative to the centroid of the convex-hull vertices.'),
    shiny::fluidRow(
      shiny::column(4,shiny::div(class='well',
        shiny::fileInput('fdiv_traits_file','1. Species traits CSV',accept='.csv'),
        shiny::fileInput('fdiv_abund_file','2. Site abundance CSV',accept='.csv'),
        shiny::uiOutput('fdiv_trait_ui'),
        shiny::actionButton('fdiv_run','Calculate 2D FDiv',class='btn-primary'),
        shiny::downloadButton('fdiv_traits_example','Example traits'),
        shiny::downloadButton('fdiv_abund_example','Example abundance'),
        shiny::p('At least three non-collinear species per site are required. Uses only the first two PCA axes.'))),
      shiny::column(8,
        shiny::h4('Validation'),shiny::verbatimTextOutput('fdiv_validation'),
        shiny::h4('Results'),shiny::tableOutput('fdiv_results'),
        shiny::downloadButton('fdiv_download','Download FDiv results CSV'),
        shiny::h4('FDiv by site'),shiny::plotOutput('fdiv_plot',height='400px'),
        shiny::h4('Method notes'),shiny::verbatimTextOutput('fdiv_notes'))))
}
ptq_fdiv_server <- function(input,output,session) {
  output$fdiv_traits_example <- shiny::downloadHandler(
    filename=function() 'fd_traits_example.csv',
    content=function(file)file.copy(file.path('data','fd_traits_example.csv'),file,overwrite=TRUE))
  output$fdiv_abund_example <- shiny::downloadHandler(
    filename=function() 'fd_abundance_example.csv',
    content=function(file)file.copy(file.path('data','fd_abundance_example.csv'),file,overwrite=TRUE))
  read_csv <- function(f) {
    shiny::req(f)
    tryCatch(utils::read.csv(f$datapath,check.names=FALSE,stringsAsFactors=FALSE,
      na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
      error=function(e)shiny::validate(shiny::need(FALSE,paste('Cannot read CSV:',conditionMessage(e)))))
  }
  dat <- shiny::reactive({
    t <- read_csv(input$fdiv_traits_file); a <- read_csv(input$fdiv_abund_file)
    names(t) <- trimws(names(t)); names(a) <- trimws(names(a))
    shiny::validate(shiny::need(!anyDuplicated(names(t)) && !anyDuplicated(names(a)),
      'Duplicate column names found.'),
      shiny::need('species_id' %in% names(t) && 'site_id' %in% names(a),
      'Traits need species_id and abundance needs site_id.'))
    t$species_id <- trimws(as.character(t$species_id))
    a$site_id <- trimws(as.character(a$site_id))
    shiny::validate(shiny::need(nrow(t)>=3 && nrow(a)>0,'At least three species and one site required.'),
      shiny::need(!anyNA(t$species_id) && all(nzchar(t$species_id)) && !anyDuplicated(t$species_id),
        'Invalid or duplicated species_id.'),
      shiny::need(!anyNA(a$site_id) && all(nzchar(a$site_id)) && !anyDuplicated(a$site_id),
        'Invalid or duplicated site_id.'))
    tc <- setdiff(names(t),'species_id'); ac <- setdiff(names(a),'site_id')
    shiny::validate(shiny::need(length(tc)>=2 && length(ac)>0,'Need two traits and species abundance columns.'))
    for(nm in tc) {
      v <- suppressWarnings(as.numeric(t[[nm]]))
      if(all(is.na(t[[nm]]) | is.finite(v))) t[[nm]] <- v
    }
    for(nm in ac) {
      v <- suppressWarnings(as.numeric(a[[nm]]))
      shiny::validate(shiny::need(all(is.na(a[[nm]]) | is.finite(v)),paste('Non-numeric abundance:',nm)))
      a[[nm]] <- v
    }
    m <- as.matrix(a[ac])
    shiny::validate(shiny::need(!anyNA(m),'Correct missing abundance cells before calculation.'),
      shiny::need(all(is.finite(m)) && all(m>=0),'Abundance must be finite and nonnegative.'))
    numeric_traits <- tc[vapply(t[tc],is.numeric,logical(1))]
    shiny::validate(shiny::need(length(numeric_traits)>=2,'At least two numeric traits required.'))
    list(t=t,a=a,ac=ac,m=m,numeric_traits=numeric_traits)
  })
  output$fdiv_trait_ui <- shiny::renderUI({
    x <- dat()
    shiny::selectInput('fdiv_traits','3. Choose numeric traits',
      choices=x$numeric_traits,selected=x$numeric_traits,multiple=TRUE)
  })
  output$fdiv_validation <- shiny::renderText({
    if(is.null(input$fdiv_traits_file)||is.null(input$fdiv_abund_file))
      return('Upload both CSV files or use the examples.')
    x <- dat(); chosen <- intersect(input$fdiv_traits,x$numeric_traits)
    present <- x$ac[colSums(x$m)>0]
    missing <- setdiff(present,x$t$species_id)
    paste0('Sites: ',nrow(x$a),' | Trait species: ',nrow(x$t),
      '\nSelected numeric traits: ',length(chosen),
      '\nSpecies with abundance but no trait record: ',
      if(length(missing))paste(missing,collapse=', ') else 'none',
      '\nMissing selected trait values: ',if(length(chosen))sum(is.na(x$t[chosen])) else 0)
  })
  calc <- shiny::eventReactive(input$fdiv_run,{
    x <- dat(); chosen <- intersect(input$fdiv_traits,x$numeric_traits)
    shiny::validate(shiny::need(length(chosen)>=2,'Choose at least two numeric traits.'))
    present <- x$ac[colSums(x$m)>0]
    missing <- setdiff(present,x$t$species_id)
    shiny::validate(shiny::need(!length(missing),paste('Missing trait records:',paste(missing,collapse=', '))),
      shiny::need(length(present)>=3,'At least three species must be present across sites.'))
    t <- x$t[match(present,x$t$species_id),chosen,drop=FALSE]
    shiny::validate(shiny::need(!anyNA(t) && all(is.finite(as.matrix(t))),
      'Correct missing or nonfinite trait values for species present in abundance.'),
      shiny::need(all(vapply(t,stats::sd,numeric(1))>0),'Remove constant traits.'))
    fit <- stats::prcomp(t,center=TRUE,scale.=TRUE)
    shiny::validate(shiny::need(ncol(fit$x)>=2 && fit$sdev[2]>1e-10,
      'At least two independent trait-space axes required.'))
    xy <- fit$x[,1:2,drop=FALSE]
    m <- x$m[,match(present,x$ac),drop=FALSE]
    vals <- vapply(seq_len(nrow(m)),function(i)ptq_fdiv_value(xy,m[i,]),numeric(1))
    tab <- data.frame(site_id=x$a$site_id,species_richness=rowSums(m>0),
      FDiv_2D=vals,check.names=FALSE)
    list(tab=tab,variance=100*fit$sdev^2/sum(fit$sdev^2))
  })
  output$fdiv_results <- shiny::renderTable(calc()$tab,striped=TRUE,digits=5)
  output$fdiv_download <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_FDiv_2D.csv',
    content=function(file)utils::write.csv(calc()$tab,file,row.names=FALSE,na=''))
  output$fdiv_plot <- shiny::renderPlot({
    z <- calc()$tab
    graphics::par(mar=c(8,4.5,2,1))
    graphics::barplot(z$FDiv_2D,names.arg=z$site_id,las=2,col='#1c8061',
      ylim=c(0,1),ylab='FDiv (2D)',main='Functional divergence by site')
  },width=850,height=500,res=96)
  output$fdiv_notes <- shiny::renderText({
    z <- calc()
    paste0('FDiv uses Euclidean distances in the first two standardized trait PCA axes.',
      '\nConvex hull vertex centroid is calculated separately for each site.',
      '\nFDiv=(mean radial distance + weighted mean deviation)/(mean radial distance + weighted absolute deviation).',
      '\nPC1 explains ',round(z$variance[1],1),'%; PC2 explains ',round(z$variance[2],1),'%.',
      '\nUndefined sites: ',sum(is.na(z$tab$FDiv_2D)),
      '\nThis is 2D-projected FDiv and not a full multidimensional estimate.')
  })
}
