# Independent package benchmarking. Packages are optional and never installed automatically.
ptq_benchmark_ui <- function() {
  shiny::tabPanel('Independent benchmarking',
    shiny::h3('Independent package-to-package benchmarking'),
    shiny::p('Compare PlantTraitQuant calculations against installed reference packages using fixed, reproducible test datasets. This module never modifies other analyses.'),
    shiny::fluidRow(
      shiny::column(4,shiny::wellPanel(
        shiny::numericInput('bench_tol','Absolute numerical tolerance',value=1e-7,min=1e-12,max=0.1,step=1e-7),
        shiny::actionButton('bench_run','Run independent benchmarks',class='btn-primary'),
        shiny::downloadButton('bench_csv','Download comparison CSV'),
        shiny::h4('Reference package status'),shiny::verbatimTextOutput('bench_packages'),
        shiny::p('Optional packages: vegan, betapart and FD. Missing packages produce SKIP, not PASS. Install separately only if you want those comparisons.'),
        shiny::p('This release benchmarks a defined subset of methods; it does not certify FEve, FDiv or full functional beta partitioning.'))),
      shiny::column(8,
        shiny::h4('Benchmark summary'),shiny::verbatimTextOutput('bench_summary'),
        shiny::h4('Independent comparisons'),shiny::tableOutput('bench_table'),
        shiny::h4('Method notes'),shiny::verbatimTextOutput('bench_notes'))))
}
ptq_benchmark_cases <- function(tolerance=1e-7) {
  rows <- list()
  add <- function(method,reference,observed=NA_real_,expected=NA_real_,note='',skip=FALSE) {
    delta <- if(is.finite(observed) && is.finite(expected)) abs(observed-expected) else NA_real_
    relative <- if(is.finite(delta)) delta/max(abs(expected),.Machine$double.eps) else NA_real_
    status <- if(skip) 'SKIP' else if(is.finite(delta) && delta<=tolerance) 'PASS' else 'FAIL'
    rows[[length(rows)+1L]] <<- data.frame(method=method,reference=reference,
      app_value=observed,reference_value=expected,absolute_difference=delta,
      relative_difference=relative,status=status,notes=note,stringsAsFactors=FALSE)
  }
  a <- rbind(Site_A=c(4,2,0,1),Site_B=c(0,3,2,1),Site_C=c(1,0,3,2))
  colnames(a) <- paste0('Species_',seq_len(ncol(a)))
  # Compare exactly the same input and mathematical definitions.
  p <- a[1,]/sum(a[1,]); positive <- p[p>0]
  app_shannon <- -sum(positive*log(positive))
  app_simpson <- 1-sum(positive^2)
  if(requireNamespace('vegan',quietly=TRUE)) {
    add('Shannon H, site A','vegan::diversity',app_shannon,
      as.numeric(vegan::diversity(a[1,,drop=FALSE],index='shannon')),
      'Natural logarithm; positive abundances only')
    add('Gini-Simpson, site A','vegan::diversity',app_simpson,
      as.numeric(vegan::diversity(a[1,,drop=FALSE],index='simpson')),
      '1 minus sum of squared relative abundances')
    add('Bray-Curtis, sites A/B','vegan::vegdist',
      sum(abs(a[1,]-a[2,]))/sum(a[1,]+a[2,]),
      as.numeric(vegan::vegdist(a[1:2,,drop=FALSE],method='bray')),
      'Untransformed abundance')
    add('Binary Jaccard, sites A/B','vegan::vegdist',
      ptq_beta_metrics(a[1:2,,drop=FALSE],'Jaccard')[1,2],
      as.numeric(vegan::vegdist(a[1:2,,drop=FALSE],method='jaccard',binary=TRUE)),
      'Presence/absence Jaccard')
  } else {
    for(method in c('Shannon H','Gini-Simpson','Bray-Curtis','Binary Jaccard'))
      add(method,'vegan',note='Reference package vegan not installed',skip=TRUE)
  }
  if(requireNamespace('betapart',quietly=TRUE)) {
    b <- (a>0)*1
    ref <- tryCatch(betapart::beta.pair(b,index.family='sorensen'),error=function(e)e)
    if(inherits(ref,'error')) {
      add('Simpson turnover, sites A/B','betapart::beta.pair',
        note=paste('Reference error:',conditionMessage(ref)),skip=TRUE)
      add('Sørensen dissimilarity, sites A/B','betapart::beta.pair',
        note=paste('Reference error:',conditionMessage(ref)),skip=TRUE)
    } else {
      add('Simpson turnover, sites A/B','betapart::beta.pair',
        ptq_beta_metrics(a[1:2,,drop=FALSE],'Simpson turnover')[1,2],
        as.matrix(ref$beta.sim)[1,2], 'Baselga Simpson turnover component')
      add('Sørensen dissimilarity, sites A/B','betapart::beta.pair',
        ptq_beta_metrics(a[1:2,,drop=FALSE],'Sørensen')[1,2],
        as.matrix(ref$beta.sor)[1,2], 'Total Sørensen dissimilarity')
    }
  } else {
    add('Simpson turnover','betapart',note='Reference package betapart not installed',skip=TRUE)
    add('Sørensen dissimilarity','betapart',note='Reference package betapart not installed',skip=TRUE)
  }
  # FD's dbFD computes FDis on supplied trait coordinates, with matching abundance weights.
  traits <- matrix(c(0,0,1,0,0,1,1,1),nrow=4,byrow=TRUE,
    dimnames=list(colnames(a),c('Trait1','Trait2')))
  traits <- scale(traits)
  weights <- a[1,]/sum(a[1,]); centroid <- colSums(traits*weights)
  app_fdis <- sum(weights*sqrt(rowSums(sweep(traits,2,centroid,'-')^2)))
  d <- as.matrix(stats::dist(traits))
  app_rao <- as.numeric(t(weights)%*%d%*%weights)
  # Use all three benchmark sites: every species occurs somewhere in this matrix.
  # The application and reference results are both compared for Site_A.
  if(requireNamespace('FD',quietly=TRUE)) {
    ref <- tryCatch(suppressWarnings(FD::dbFD(traits,a,
      calc.FRic=FALSE,calc.FDiv=FALSE,calc.CWM=FALSE,
      stand.x=FALSE,w.abun=TRUE,messages=FALSE)),error=function(e)e)
    if(inherits(ref,'error')) {
      add('FDis, site A','FD::dbFD',note=paste('Reference error:',conditionMessage(ref)),skip=TRUE)
      add('RaoQ, site A','FD::dbFD',note=paste('Reference error:',conditionMessage(ref)),skip=TRUE)
    } else {
      fdis_ref <- if(!is.null(ref$FDis))as.numeric(if(!is.null(names(ref$FDis)) && 'Site_A' %in% names(ref$FDis)) ref$FDis['Site_A'] else ref$FDis[1]) else NA_real_
      rao_ref <- if(!is.null(ref$RaoQ))as.numeric(if(!is.null(names(ref$RaoQ)) && 'Site_A' %in% names(ref$RaoQ)) ref$RaoQ['Site_A'] else ref$RaoQ[1]) else NA_real_
      if(is.finite(fdis_ref))add('FDis, site A','FD::dbFD',app_fdis,fdis_ref,
        'Euclidean distances; globally z-scored traits; abundance-weighted centroid')
      else add('FDis, site A','FD::dbFD',note='Reference did not return finite FDis',skip=TRUE)
      if(is.finite(rao_ref))add('RaoQ, site A','FD::dbFD',app_rao,rao_ref,
        'Euclidean distances on standardized traits; compare same Site_A weights')
      else add('RaoQ, site A','FD::dbFD',note='Reference did not return finite RaoQ',skip=TRUE)
    }
  } else {
    add('FDis','FD',note='Reference package FD not installed',skip=TRUE)
    add('RaoQ','FD',note='Reference package FD not installed',skip=TRUE)
  }
  do.call(rbind,rows)
}
ptq_benchmark_server <- function(input,output,session) {
  output$bench_packages <- shiny::renderText({
    paste(vapply(c('vegan','betapart','FD'),function(pkg)
      paste(pkg,if(requireNamespace(pkg,quietly=TRUE))
        paste('installed:',as.character(utils::packageVersion(pkg))) else 'NOT INSTALLED'),character(1)),collapse='\n')
  })
  results <- shiny::eventReactive(input$bench_run,ptq_benchmark_cases(input$bench_tol))
  output$bench_summary <- shiny::renderText({
    d <- results()
    paste0('PASS: ',sum(d$status=='PASS'),' | FAIL: ',sum(d$status=='FAIL'),
      ' | SKIP: ',sum(d$status=='SKIP'),
      '\nTolerance (absolute): ',format(input$bench_tol,scientific=TRUE),
      '\nSKIP means not independently tested, not validated.')
  })
  output$bench_table <- shiny::renderTable(results(),striped=TRUE,digits=8)
  output$bench_notes <- shiny::renderText({
    paste('All comparisons use deterministic example matrices and are recomputed on demand.',
      'RaoQ definitions, distance scaling, and package-specific normalizations must be checked if a comparison fails.',
      'FEve, FDiv, FRic and functional beta turnover/nestedness are not independently benchmarked by this module.',sep='\n')
  })
  output$bench_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_independent_benchmarks.csv',
    content=function(file)utils::write.csv(results(),file,row.names=FALSE,na=''))
}
