# Pairwise expected between-community trait dissimilarity.
# Explicitly descriptive: not a turnover/nestedness decomposition.
ptq_pairwise_ui <- function() {
  controls <- shiny::column(4,
    shiny::wellPanel(
      shiny::fileInput('pair_traits_file','Species traits CSV',accept='.csv'),
      shiny::fileInput('pair_abundance_file','Site abundance CSV',accept='.csv'),
      shiny::uiOutput('pair_trait_selector'),
      shiny::checkboxInput('pair_binary','Use presence/absence instead of abundance',FALSE),
      shiny::actionButton('pair_calculate','Calculate pairwise dissimilarity',class='btn-primary'),
      shiny::downloadButton('pair_example_traits','Example traits'),
      shiny::downloadButton('pair_example_abundance','Example abundance')))
  results <- shiny::column(8,
    shiny::h4('Data checks'),shiny::verbatimTextOutput('pair_checks'),
    shiny::h4('Between-community expected trait dissimilarity'),
    shiny::tableOutput('pair_matrix'),
    shiny::downloadButton('pair_matrix_download','Download matrix CSV'),
    shiny::h4('Site-level expected within-community dissimilarity'),
    shiny::tableOutput('pair_within'),
    shiny::downloadButton('pair_within_download','Download within-site CSV'),
    shiny::h4('Method notes'),shiny::verbatimTextOutput('pair_notes'))
  shiny::tabPanel('Pairwise functional differences',
    shiny::h3('Abundance-weighted pairwise functional dissimilarity'),
    shiny::p('Calculate expected trait dissimilarity between two independently drawn species, one from each site. This is a descriptive cross-community quantity, not a formal turnover/nestedness partition.'),
    shiny::fluidRow(controls,results))
}
ptq_pairwise_server <- function(input,output,session) {
  output$pair_example_traits <- shiny::downloadHandler(
    filename=function() 'fd_traits_example.csv',
    content=function(file) file.copy(file.path('data','fd_traits_example.csv'),file,overwrite=TRUE))
  output$pair_example_abundance <- shiny::downloadHandler(
    filename=function() 'fd_abundance_example.csv',
    content=function(file) file.copy(file.path('data','fd_abundance_example.csv'),file,overwrite=TRUE))
  read_csv <- function(file) {
    shiny::req(file)
    tryCatch(utils::read.csv(file$datapath,check.names=FALSE,stringsAsFactors=FALSE,
       na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
       error=function(e) shiny::validate(shiny::need(FALSE,conditionMessage(e))))
  }
  data_pair <- shiny::reactive({
    t <- read_csv(input$pair_traits_file)
    a <- read_csv(input$pair_abundance_file)
    names(t) <- trimws(names(t));names(a) <- trimws(names(a))
    shiny::validate(
      shiny::need(!anyDuplicated(names(t)) && !anyDuplicated(names(a)),'Duplicate column names.'),
      shiny::need('species_id' %in% names(t),'Trait file needs species_id.'),
      shiny::need('site_id' %in% names(a),'Abundance file needs site_id.'),
      shiny::need(nrow(a)>=2 && nrow(a)<=200,'Provide between 2 and 200 sites.'))
    t$species_id <- trimws(as.character(t$species_id))
    a$site_id <- trimws(as.character(a$site_id))
    shiny::validate(
      shiny::need(all(!is.na(t$species_id) & nzchar(t$species_id)) && !anyDuplicated(t$species_id),'Invalid or repeated species_id.'),
      shiny::need(all(!is.na(a$site_id) & nzchar(a$site_id)) && !anyDuplicated(a$site_id),'Invalid or repeated site_id.'))
    species <- setdiff(names(a),'site_id')
    shiny::validate(shiny::need(length(species)>0,'Abundance needs species columns.'))
    m <- sapply(a[species],function(x) suppressWarnings(as.numeric(x)))
    m <- matrix(m,nrow=nrow(a),ncol=length(species),dimnames=list(a$site_id,species))
    shiny::validate(shiny::need(all(is.finite(m) & m>=0),'All abundance cells must be finite, nonnegative numeric values.'))
    shiny::validate(shiny::need(all(rowSums(m)>0),'Empty sites are not allowed.'))
    tc <- setdiff(names(t),'species_id')
    for(nm in tc) {
      raw <- t[[nm]]; converted <- suppressWarnings(as.numeric(raw))
      if(all(is.na(raw) | is.finite(converted))) t[[nm]] <- converted
    }
    numeric <- tc[vapply(t[tc],is.numeric,logical(1))]
    shiny::validate(shiny::need(length(numeric)>=1,'At least one numeric trait required.'))
    list(t=t,a=a,m=m,numeric=numeric)
  })
  output$pair_trait_selector <- shiny::renderUI({
    x <- data_pair()
    shiny::selectInput('pair_traits','Select numerical traits',choices=x$numeric,
      selected=x$numeric,multiple=TRUE)
  })
  output$pair_checks <- shiny::renderText({
    if(is.null(input$pair_traits_file) || is.null(input$pair_abundance_file))
      return('Upload both files or use the example pair.')
    x <- data_pair()
    used <- colnames(x$m)[colSums(x$m)>0]
    unmatched <- setdiff(used,x$t$species_id)
    paste('Sites:',nrow(x$a),'| Species observed:',length(used),
      '| Numeric traits:',length(x$numeric),
      '\nObserved species missing trait records:',
      if(length(unmatched))paste(unmatched,collapse=', ') else 'none')
  })
  calculated <- shiny::eventReactive(input$pair_calculate,{
    x <- data_pair()
    chosen <- intersect(input$pair_traits,x$numeric)
    shiny::validate(shiny::need(length(chosen)>=1,'Choose at least one numeric trait.'))
    used <- colnames(x$m)[colSums(x$m)>0]
    absent <- setdiff(used,x$t$species_id)
    shiny::validate(shiny::need(length(absent)==0,
      paste('Missing traits for species:',paste(absent,collapse=', '))))
    t <- x$t[match(used,x$t$species_id),chosen,drop=FALSE]
    shiny::validate(shiny::need(all(stats::complete.cases(t)) && all(is.finite(as.matrix(t))),
      'Correct missing or infinite selected traits for observed species.'))
    sdv <- vapply(t,stats::sd,numeric(1))
    shiny::validate(shiny::need(all(is.finite(sdv) & sdv>0),
      'Remove constant traits; selected traits need variation.'))
    scaled <- scale(as.matrix(t))
    d <- as.matrix(stats::dist(scaled))
    max_d <- max(d)
    shiny::validate(shiny::need(is.finite(max_d) && max_d>0,
      'No variation in species trait distances.'))
    d <- d/max_d
    a <- x$m[,used,drop=FALSE]
    if(isTRUE(input$pair_binary)) a <- (a>0)*1
    p <- a/rowSums(a)
    # E[d(X,Y)] where X and Y are independently sampled from two communities.
    # Diagonal is within-site expected dissimilarity, not forced to zero.
    between <- p %*% d %*% t(p)
    dimnames(between) <- list(rownames(a),rownames(a))
    within <- data.frame(site_id=rownames(a),
      within_expected_dissimilarity=diag(between),row.names=NULL)
    list(between=between,within=within)
  })
  output$pair_matrix <- shiny::renderTable({
    z <- calculated()$between
    as.data.frame(round(z[seq_len(min(nrow(z),12)),seq_len(min(ncol(z),12)),drop=FALSE],4))
  },rownames=TRUE,striped=TRUE)
  output$pair_within <- shiny::renderTable(calculated()$within,striped=TRUE)
  output$pair_matrix_download <- shiny::downloadHandler(
    filename=function() 'PlantTraitQuant_cross_community_trait_dissimilarity.csv',
    content=function(file) utils::write.csv(calculated()$between,file,row.names=TRUE))
  output$pair_within_download <- shiny::downloadHandler(
    filename=function() 'PlantTraitQuant_within_community_trait_dissimilarity.csv',
    content=function(file) utils::write.csv(calculated()$within,file,row.names=FALSE))
  output$pair_notes <- shiny::renderText({
    shiny::req(calculated())
    paste('Selected traits are z-standardized across observed species.',
      'Euclidean species distances are divided by their maximum.',
      'Each matrix cell is sum_i sum_j p(siteA,i) p(siteB,j) d(i,j).',
      'The diagonal gives within-community expected dissimilarity and may be positive.',
      'This matrix is NOT a distance metric: do not apply PCoA or clustering to it as though its diagonal were zero.',
      'It is not a functional turnover/nestedness decomposition.',sep='\n')
  })
}
