# Consolidated functional diversity dashboard, base R + shiny only.
ptq_dashboard_ui <- function() {
  shiny::tabPanel('FD dashboard (all indices)',
    shiny::h3('Consolidated functional diversity dashboard'),
    shiny::p('One pair of files and one trait selection for FDis, RaoQ, 2D FRic, FEve and 2D FDiv.'),
    shiny::fluidRow(
      shiny::column(4, shiny::div(class='well',
        shiny::fileInput('dash_traits_file','1. Species traits CSV',accept='.csv'),
        shiny::fileInput('dash_abund_file','2. Site abundance CSV',accept='.csv'),
        shiny::uiOutput('dash_trait_selector'),
        shiny::actionButton('dash_run','Calculate all indices',class='btn-primary'),
        shiny::downloadButton('dash_example_traits','Example species traits'),
        shiny::downloadButton('dash_example_abund','Example site abundance'),
        shiny::p('Numeric traits only. All indices use the same standardized trait space. FRic and FDiv use the first two PCA axes.'))),
      shiny::column(8,
        shiny::h4('Data quality and matching'),shiny::verbatimTextOutput('dash_validation'),
        shiny::h4('All site-level indices'),shiny::tableOutput('dash_results'),
        shiny::downloadButton('dash_download','Download all indices CSV'),
        shiny::h4('Compare indices across sites'),
        shiny::selectInput('dash_index','Index to display',
          choices=c('FDis','RaoQ','FRic_2D','FEve','FDiv_2D')),
        shiny::plotOutput('dash_plot',height='400px'),
        shiny::h4('Calculation notes'),shiny::verbatimTextOutput('dash_notes'),
        shiny::hr(),
        shiny::h3('Visual exploration of functional diversity'),
        shiny::p('These plots use the calculated dashboard results. No additional data upload or packages are needed.'),
        shiny::h4('Standardized site-by-index heatmap'),
        shiny::plotOutput('dash_heatmap',height='460px'),
        shiny::downloadButton('dash_heatmap_png','Download heatmap PNG'),
        shiny::h4('Correlation among indices'),
        shiny::verbatimTextOutput('dash_corr_note'),
        shiny::plotOutput('dash_corr_plot',height='430px'),
        shiny::downloadButton('dash_corr_csv','Download correlations CSV'),
        shiny::h4('Hierarchical site clustering'),
        shiny::plotOutput('dash_cluster_plot',height='430px'),
        shiny::downloadButton('dash_cluster_csv','Download clustering input CSV'))))
}
ptq_dashboard_server <- function(input,output,session) {
  output$dash_example_traits <- shiny::downloadHandler(
    filename=function() 'fd_traits_example.csv',
    content=function(file) file.copy(file.path('data','fd_traits_example.csv'),file,overwrite=TRUE))
  output$dash_example_abund <- shiny::downloadHandler(
    filename=function() 'fd_abundance_example.csv',
    content=function(file) file.copy(file.path('data','fd_abundance_example.csv'),file,overwrite=TRUE))
  read_csv <- function(upload) {
    shiny::req(upload)
    tryCatch(utils::read.csv(upload$datapath,check.names=FALSE,
      stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
      error=function(e) shiny::validate(shiny::need(FALSE,paste('CSV read error:',conditionMessage(e)))))
  }
  datasets <- shiny::reactive({
    t <- read_csv(input$dash_traits_file)
    a <- read_csv(input$dash_abund_file)
    names(t) <- trimws(names(t)); names(a) <- trimws(names(a))
    shiny::validate(
      shiny::need(!anyDuplicated(names(t)) && !anyDuplicated(names(a)),'Duplicate column names found.'),
      shiny::need('species_id' %in% names(t),'Trait file must contain species_id.'),
      shiny::need('site_id' %in% names(a),'Abundance file must contain site_id.'),
      shiny::need(nrow(t)>=3 && nrow(a)>0,'Need at least 3 trait species and 1 site.'))
    t$species_id <- trimws(as.character(t$species_id))
    a$site_id <- trimws(as.character(a$site_id))
    shiny::validate(
      shiny::need(!anyNA(t$species_id) && all(nzchar(t$species_id)) && !anyDuplicated(t$species_id),'Invalid or duplicate species_id.'),
      shiny::need(!anyNA(a$site_id) && all(nzchar(a$site_id)) && !anyDuplicated(a$site_id),'Invalid or duplicate site_id.'))
    tc <- setdiff(names(t),'species_id')
    ac <- setdiff(names(a),'site_id')
    shiny::validate(shiny::need(length(tc)>=2 && length(ac)>0,'Need at least 2 traits and one species abundance column.'))
    for(nm in tc) {
      v <- suppressWarnings(as.numeric(t[[nm]]))
      if(all(is.na(t[[nm]]) | is.finite(v))) t[[nm]] <- v
    }
    for(nm in ac) {
      v <- suppressWarnings(as.numeric(a[[nm]]))
      shiny::validate(shiny::need(all(is.na(a[[nm]]) | is.finite(v)),paste('Non-numeric abundance in:',nm)))
      a[[nm]] <- v
    }
    m <- as.matrix(a[ac])
    shiny::validate(shiny::need(!anyNA(m),'Missing abundance values: correct the file.'),
      shiny::need(all(is.finite(m) & m>=0),'Abundance values must be finite and nonnegative.'))
    nums <- tc[vapply(t[tc],is.numeric,logical(1))]
    shiny::validate(shiny::need(length(nums)>=2,'At least 2 numeric traits are required.'))
    list(t=t,a=a,m=m,ac=ac,nums=nums)
  })
  output$dash_trait_selector <- shiny::renderUI({
    x <- datasets()
    shiny::selectInput('dash_traits','3. Select numeric traits',choices=x$nums,
      selected=x$nums,multiple=TRUE)
  })
  output$dash_validation <- shiny::renderText({
    if(is.null(input$dash_traits_file) || is.null(input$dash_abund_file))
      return('Upload both files or download the paired examples.')
    x <- datasets()
    present <- x$ac[colSums(x$m)>0]
    missing <- setdiff(present,x$t$species_id)
    selected <- intersect(input$dash_traits,x$nums)
    paste0('Sites: ',nrow(x$a),' | Trait species: ',nrow(x$t),
      '\nSpecies occurring in sites: ',length(present),
      '\nSelected numeric traits: ',length(selected),
      '\nPresent species without trait records: ',
      if(length(missing)) paste(missing,collapse=', ') else 'none',
      '\nMissing selected trait values in present species: ',
      if(length(selected) && length(present)) {
        matched <- match(intersect(present,x$t$species_id),x$t$species_id)
        sum(is.na(x$t[matched,selected,drop=FALSE]))
      } else 0,
      '\nSites with zero total abundance: ',sum(rowSums(x$m)==0))
  })
  calculated <- shiny::eventReactive(input$dash_run,{
    x <- datasets()
    selected <- intersect(input$dash_traits,x$nums)
    shiny::validate(shiny::need(length(selected)>=2,'Select at least two numeric traits.'))
    present <- x$ac[colSums(x$m)>0]
    missing <- setdiff(present,x$t$species_id)
    shiny::validate(shiny::need(length(missing)==0,
      paste('Species without trait records:',paste(missing,collapse=', '))),
      shiny::need(length(present)>=3,'At least three species must occur across all sites.'))
    t <- x$t[match(present,x$t$species_id),selected,drop=FALSE]
    shiny::validate(shiny::need(!anyNA(t) && all(is.finite(as.matrix(t))),
      'Missing or nonfinite trait values for species occurring in sites.'))
    variability <- vapply(t,stats::sd,numeric(1))
    shiny::validate(shiny::need(all(is.finite(variability) & variability>0),
      'Remove constant selected traits.'))
    z <- scale(as.matrix(t),center=TRUE,scale=TRUE)
    dm <- as.matrix(stats::dist(z))
    fit <- stats::prcomp(z,center=FALSE,scale.=FALSE)
    shiny::validate(shiny::need(ncol(fit$x)>=2 && fit$sdev[2]>1e-10,
      'Two independent PCA axes are required for FRic_2D and FDiv_2D.'))
    xy <- fit$x[,1:2,drop=FALSE]
    m <- x$m[,match(present,x$ac),drop=FALSE]
    site_metrics <- lapply(seq_len(nrow(m)),function(i) {
      abundance <- as.numeric(m[i,])
      richness <- sum(abundance>0)
      if(sum(abundance)==0) return(c(richness=0,FDis=NA_real_,RaoQ=NA_real_,
                                      FRic_2D=NA_real_,FEve=NA_real_,FDiv_2D=NA_real_))
      p <- abundance/sum(abundance)
      centroid <- colSums(z*p)
      radial <- sqrt(rowSums(sweep(z,2,centroid,'-')^2))
      c(richness=richness,
        FDis=sum(p*radial),
        RaoQ=as.numeric(crossprod(p,dm%*%p)),
        FRic_2D=ptq_fric2d(xy,abundance),
        FEve=ptq_feve_value(dm,abundance),
        FDiv_2D=ptq_fdiv_value(xy,abundance))
    })
    results <- data.frame(site_id=x$a$site_id,do.call(rbind,site_metrics),
                          row.names=NULL,check.names=FALSE)
    list(results=results,variance=100*fit$sdev^2/sum(fit$sdev^2),
         species=length(present),traits=selected)
  })
  output$dash_results <- shiny::renderTable(calculated()$results,striped=TRUE,digits=5)
  output$dash_download <- shiny::downloadHandler(
    filename=function() 'PlantTraitQuant_All_Functional_Diversity.csv',
    content=function(file) utils::write.csv(calculated()$results,file,row.names=FALSE,na=''))
  output$dash_plot <- shiny::renderPlot({
    z <- calculated()$results
    index <- input$dash_index
    shiny::req(index)
    y <- z[[index]]
    shiny::validate(shiny::need(any(is.finite(y)),paste('No defined',index,'values to plot.')))
    graphics::par(mar=c(8,4.5,2,1))
    graphics::barplot(y,names.arg=z$site_id,las=2,col='#237d60',
      ylab=index,main=paste(index,'across sites'))
  },width=850,height=480,res=96)
  output$dash_notes <- shiny::renderText({
    z <- calculated()
    paste0('Common standardized trait matrix: ',z$species,' species, ',
      length(z$traits),' selected traits.',
      '\nFDis: abundance-weighted mean Euclidean distance to weighted centroid.',
      '\nRaoQ: sum of p_i p_j times Euclidean trait distance.',
      '\nFRic_2D: convex-hull area of present species on first two PCA axes.',
      '\nFEve: minimum spanning tree based abundance-weighted evenness.',
      '\nFDiv_2D: abundance-weighted divergence relative to hull-vertex centroid.',
      '\nPC1 explained variance: ',round(z$variance[1],2),'%; PC2: ',round(z$variance[2],2),'%.',
      '\nNA means index undefined for that site, often because richness or geometry is insufficient.',
      '\n2D indices are not full-dimensional equivalents. Index scales differ; do not directly compare raw magnitudes.',
      '\nThese implementations require independent numerical benchmarking before publication.')
  })

  # All visualization methods use dashboard output, with no new dependencies.
  metric_names <- c('FDis','RaoQ','FRic_2D','FEve','FDiv_2D')
  visual_data <- shiny::reactive({
    tab <- calculated()$results
    m <- as.matrix(tab[,metric_names,drop=FALSE])
    storage.mode(m) <- 'double'
    rownames(m) <- tab$site_id
    usable <- vapply(seq_len(ncol(m)),function(j) {
      v <- m[,j]
      sum(is.finite(v))>=2 && stats::sd(v[is.finite(v)])>0
    },logical(1))
    m <- m[,usable,drop=FALSE]
    shiny::validate(shiny::need(ncol(m)>0,
      'At least one index must have two or more distinct finite site values.'))
    z <- apply(m,2,function(v) {
      mu <- mean(v[is.finite(v)])
      sig <- stats::sd(v[is.finite(v)])
      (v-mu)/sig
    })
    z <- matrix(z,nrow=nrow(m),ncol=ncol(m),
      dimnames=dimnames(m))
    list(raw=m,z=z)
  })
  draw_heatmap <- function() {
    v <- visual_data()
    m <- v$z
    # Explicitly mark undefined values as blank rather than imputing.
    nr <- nrow(m); nc <- ncol(m)
    graphics::par(mar=c(8,10,3,2))
    palette <- grDevices::colorRampPalette(c('#254a85','#e9edf1','#a54a35'))(101)
    bounds <- max(1,max(abs(m),na.rm=TRUE))
    plot(NA,xlim=c(.5,nc+.5),ylim=c(nr+.5,.5),
      xlab='',ylab='',xaxt='n',yaxt='n',bty='n',
      main='Functional diversity heatmap (within-index z-scores)')
    for(i in seq_len(nr)) for(j in seq_len(nc)) {
      v0 <- m[i,j]
      fill <- if(is.finite(v0)) {
        k <- max(1,min(101,round((v0/bounds+1)*50)+1))
        palette[k]
      } else '#d9d9d9'
      graphics::rect(j-.48,i-.48,j+.48,i+.48,col=fill,border='white')
    }
    graphics::axis(1,at=seq_len(nc),labels=colnames(m),las=2,cex.axis=.85)
    graphics::axis(2,at=seq_len(nr),labels=rownames(m),las=2,cex.axis=.75)
    graphics::mtext('Grey = undefined; each index standardized separately',
      side=1,line=6,cex=.7)
  }
  output$dash_heatmap <- shiny::renderPlot(draw_heatmap(),
    width=900,height=560,res=96)
  output$dash_heatmap_png <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_FD_heatmap.png',
    content=function(file) {
      grDevices::png(file,width=1100,height=700,res=110)
      on.exit(grDevices::dev.off(),add=TRUE)
      draw_heatmap()
    })
  correlation <- shiny::reactive({
    m <- visual_data()$raw
    if(ncol(m)<2) return(NULL)
    stats::cor(m,use='pairwise.complete.obs',method='pearson')
  })
  output$dash_corr_note <- shiny::renderText({
    m <- visual_data()$raw
    paste('Pearson correlations among',ncol(m),'variable indices across',
      nrow(m),'sites. Each pair uses sites with both values observed.',
      'At least two variable indices are required. Correlation does not imply causation.')
  })
  output$dash_corr_plot <- shiny::renderPlot({
    r <- correlation()
    shiny::validate(shiny::need(!is.null(r),'At least two variable indices are required.'))
    n <- ncol(r)
    graphics::par(mar=c(8,9,2,2))
    pal <- grDevices::colorRampPalette(c('#315b9a','#ffffff','#b34436'))(101)
    plot(NA,xlim=c(.5,n+.5),ylim=c(n+.5,.5),
      xaxt='n',yaxt='n',xlab='',ylab='',bty='n',
      main='Pearson correlation between FD indices')
    for(i in seq_len(n)) for(j in seq_len(n)) {
      val <- r[i,j]
      k <- if(is.finite(val))max(1,min(101,round((val+1)*50)+1)) else NA_integer_
      graphics::rect(j-.48,i-.48,j+.48,i+.48,
        col=if(is.na(k))'#d9d9d9' else pal[k],border='white')
      graphics::text(j,i,if(is.finite(val))sprintf('%.2f',val) else 'NA',cex=.8)
    }
    graphics::axis(1,at=seq_len(n),labels=colnames(r),las=2)
    graphics::axis(2,at=seq_len(n),labels=rownames(r),las=2)
  },width=800,height=550,res=96)
  output$dash_corr_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_FD_correlations.csv',
    content=function(file) {
      r <- correlation()
      shiny::validate(shiny::need(!is.null(r),'Need two variable indices.'))
      utils::write.csv(r,file,row.names=TRUE,na='')
    })
  cluster_input <- shiny::reactive({
    z <- visual_data()$z
    # Complete cases only: no automatic imputation of undefined FD indices.
    z <- z[stats::complete.cases(z),,drop=FALSE]
    shiny::validate(shiny::need(nrow(z)>=2,
      'Clustering requires at least two sites with all selected variable indices defined.'))
    z
  })
  output$dash_cluster_plot <- shiny::renderPlot({
    z <- cluster_input()
    hc <- stats::hclust(stats::dist(z),method='average')
    graphics::par(mar=c(4.5,4.5,3,2))
    graphics::plot(hc,main='Site clustering (average linkage)',
      xlab='Sites',sub='Euclidean distance on standardized FD indices',cex=.8)
  },width=900,height=540,res=96)
  output$dash_cluster_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_FD_clustering_input.csv',
    content=function(file)utils::write.csv(cluster_input(),file,row.names=TRUE))

}
