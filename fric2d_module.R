
# Two-dimensional convex-hull functional richness, a deliberately bounded implementation.
ptq_fric2d <- function(coords, abund) {
  idx <- which(abund > 0)
  if (length(idx) < 3L) return(NA_real_)
  x <- coords[idx, , drop=FALSE]
  if (qr(cbind(1,x))$rank < 3L) return(NA_real_)
  hull <- grDevices::chull(x[,1], x[,2])
  polygon <- x[c(hull,hull[1]), ,drop=FALSE]
  abs(sum(polygon[-nrow(polygon),1]*polygon[-1,2] -
          polygon[-1,1]*polygon[-nrow(polygon),2]))/2
}
ptq_fric_ui <- function() {
  shiny::tabPanel('Functional richness (2D FRic)',
    shiny::h3('Functional richness: two-dimensional convex-hull area'),
    shiny::p('A deliberately limited, transparent FRic implementation for testing. Upload species traits and site abundance independently of other tabs.'),
    shiny::fluidRow(
      shiny::column(4,shiny::div(class='well',
        shiny::fileInput('fric_traits_file','1. Species traits CSV',accept='.csv'),
        shiny::fileInput('fric_abund_file','2. Site abundance CSV',accept='.csv'),
        shiny::uiOutput('fric_trait_ui'),
        shiny::actionButton('fric_run','Calculate 2D FRic',class='btn-primary'),
        shiny::downloadButton('fric_traits_example','Example traits'),
        shiny::downloadButton('fric_abund_example','Example abundance'),
        shiny::p('The selected numeric traits are standardized globally, then reduced to two PCA axes. Convex-hull area is computed per site. This is 2D FRic, not full multidimensional FRic.'))),
      shiny::column(8,
        shiny::h4('Input validation'),shiny::verbatimTextOutput('fric_validation'),
        shiny::h4('Site-level 2D FRic'),shiny::tableOutput('fric_table'),
        shiny::downloadButton('fric_csv','Download FRic CSV'),
        shiny::h4('Trait-space ordination'),shiny::plotOutput('fric_plot',height='420px'),
        shiny::h4('Calculation notes'),shiny::verbatimTextOutput('fric_notes'))))
}
ptq_fric_server <- function(input,output,session) {
  output$fric_traits_example <- shiny::downloadHandler(
    filename=function()'fd_traits_example.csv',
    content=function(file)file.copy(file.path('data','fd_traits_example.csv'),file,overwrite=TRUE))
  output$fric_abund_example <- shiny::downloadHandler(
    filename=function()'fd_abundance_example.csv',
    content=function(file)file.copy(file.path('data','fd_abundance_example.csv'),file,overwrite=TRUE))
  read_input <- function(f) {
    shiny::req(f)
    tryCatch(utils::read.csv(f$datapath,check.names=FALSE,stringsAsFactors=FALSE,
               na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
      error=function(e) shiny::validate(shiny::need(FALSE,paste('Cannot read CSV:',conditionMessage(e)))))
  }
  inputs <- shiny::reactive({
    t <- read_input(input$fric_traits_file)
    a <- read_input(input$fric_abund_file)
    names(t) <- trimws(names(t)); names(a) <- trimws(names(a))
    shiny::validate(shiny::need(!anyDuplicated(names(t)) && !anyDuplicated(names(a)),
      'Duplicate column names found.'),
      shiny::need('species_id' %in% names(t) && 'site_id' %in% names(a),
      'Traits need species_id and abundance needs site_id.'))
    t$species_id <- trimws(as.character(t$species_id))
    a$site_id <- trimws(as.character(a$site_id))
    shiny::validate(shiny::need(nrow(t)>=3 && nrow(a)>0,'At least three species and one site required.'),
      shiny::need(!anyNA(t$species_id) && all(nzchar(t$species_id)) &&
       !anyDuplicated(t$species_id),'Invalid or duplicated species_id.'),
      shiny::need(!anyNA(a$site_id) && all(nzchar(a$site_id)) &&
       !anyDuplicated(a$site_id),'Invalid or duplicated site_id.'))
    tc <- setdiff(names(t),'species_id')
    ac <- setdiff(names(a),'site_id')
    shiny::validate(shiny::need(length(tc)>=2 && length(ac)>0,'Need two traits and abundance columns.'))
    for(nm in tc) {
      v <- suppressWarnings(as.numeric(t[[nm]]))
      if(all(is.na(t[[nm]]) | is.finite(v))) t[[nm]] <- v
    }
    for(nm in ac) {
      v <- suppressWarnings(as.numeric(a[[nm]]))
      shiny::validate(shiny::need(all(is.na(a[[nm]]) | is.finite(v)),
        paste('Non-numeric abundance in',nm)))
      a[[nm]] <- v
    }
    m <- as.matrix(a[ac])
    shiny::validate(shiny::need(!anyNA(m),'Missing abundance cells must be corrected.'),
      shiny::need(all(is.finite(m)) && all(m>=0),'Abundances must be finite and nonnegative.'))
    numeric_traits <- tc[vapply(t[tc],is.numeric,logical(1))]
    shiny::validate(shiny::need(length(numeric_traits)>=2,'At least two numeric traits required.'))
    list(t=t,a=a,ac=ac,numeric_traits=numeric_traits,m=m)
  })
  output$fric_trait_ui <- shiny::renderUI({
    x <- inputs()
    shiny::selectInput('fric_traits','3. Choose numeric traits',
      choices=x$numeric_traits,selected=x$numeric_traits,multiple=TRUE)
  })
  output$fric_validation <- shiny::renderText({
    if(is.null(input$fric_traits_file)||is.null(input$fric_abund_file))
      return('Upload two files or use the examples.')
    x <- inputs()
    present <- x$ac[colSums(x$m)>0]
    unmatched <- setdiff(present,x$t$species_id)
    chosen <- intersect(input$fric_traits,x$numeric_traits)
    paste0('Sites: ',nrow(x$a),' | Trait species: ',nrow(x$t),
      '\nSelected numeric traits: ',length(chosen),
      '\nPositive-abundance species without traits: ',
      if(length(unmatched))paste(unmatched,collapse=', ') else 'none',
      '\nMissing selected trait cells among trait records: ',
      if(length(chosen))sum(is.na(x$t[chosen])) else 0,
      '\nFRic is undefined for fewer than 3 distinct non-collinear species in 2D.')
  })
  calc <- shiny::eventReactive(input$fric_run,{
    x <- inputs()
    chosen <- intersect(input$fric_traits,x$numeric_traits)
    shiny::validate(shiny::need(length(chosen)>=2,'Choose at least two numeric traits.'))
    present <- x$ac[colSums(x$m)>0]
    missing_species <- setdiff(present,x$t$species_id)
    shiny::validate(shiny::need(length(missing_species)==0,
      paste('Missing trait records for:',paste(missing_species,collapse=', '))),
      shiny::need(length(present)>=3,'At least three species must occur across sites.'))
    t <- x$t[match(present,x$t$species_id),chosen,drop=FALSE]
    shiny::validate(shiny::need(!anyNA(t) && all(is.finite(as.matrix(t))),
      'Missing or nonfinite selected trait values for present species.'))
    sdv <- vapply(t,stats::sd,numeric(1))
    shiny::validate(shiny::need(all(is.finite(sdv) & sdv>0),
      'Remove constant selected traits.'))
    z <- scale(as.matrix(t))
    fit <- stats::prcomp(z,center=FALSE,scale.=FALSE)
    shiny::validate(shiny::need(ncol(fit$x)>=2 && fit$sdev[2]>1e-10,
      'The trait matrix has fewer than two independent axes.'))
    xy <- fit$x[,1:2,drop=FALSE]
    rownames(xy) <- present
    m <- x$m[,match(present,x$ac),drop=FALSE]
    values <- vapply(seq_len(nrow(m)),function(i)ptq_fric2d(xy,m[i,]),numeric(1))
    result <- data.frame(site_id=x$a$site_id,
      species_richness=rowSums(m>0),
      FRic_2D=values,check.names=FALSE)
    list(table=result,xy=xy,variance=100*fit$sdev^2/sum(fit$sdev^2))
  })
  output$fric_table <- shiny::renderTable(calc()$table,striped=TRUE,digits=5)
  output$fric_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_FRic_2D.csv',
    content=function(file)utils::write.csv(calc()$table,file,row.names=FALSE,na=''))
  output$fric_plot <- shiny::renderPlot({
    z <- calc();xy <- z$xy
    graphics::par(mar=c(4.5,4.5,2.5,1.5))
    graphics::plot(xy[,1],xy[,2],pch=19,col='#176d51',
      xlab=paste0('PC1 (',round(z$variance[1],1),'%)'),
      ylab=paste0('PC2 (',round(z$variance[2],1),'%)'),
      main='Shared standardized species trait space')
    if(nrow(xy)<=60) graphics::text(xy[,1],xy[,2],labels=rownames(xy),pos=3,cex=.7)
    graphics::abline(h=0,v=0,lty=3,col='grey70')
  },width=850,height=530,res=96)
  output$fric_notes <- shiny::renderText({
    z <- calc()
    paste0('FRic is area of the convex hull in the first two axes of globally standardized trait PCA.',
      '\nPC1 variance: ',round(z$variance[1],2),'% | PC2 variance: ',round(z$variance[2],2),'%.',
      '\nUndefined sites (fewer than 3 non-collinear present species): ',
      sum(is.na(z$table$FRic_2D)),
      '\nFRic_2D ignores species abundance magnitudes, using presence only.',
      '\nArea units are squared standardized PCA-score units.',
      '\nDo not compare these values directly with full-dimensional FRic.')
  })
}
