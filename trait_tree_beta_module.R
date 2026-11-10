# Trait-dendrogram functional beta partition (incidence based).
# Branch lengths are calculated from average-linkage clustering of standardized
# Euclidean species trait distances. These are dendrogram branch lengths, not
# multidimensional convex-hull functional beta measures.
ptq_tree_branch_matrix <- function(hc, species) {
  n <- length(species)
  members <- vector('list', n-1L)
  lengths <- numeric(2L*n-2L)
  incidence <- matrix(FALSE,nrow=2L*n-2L,ncol=n)
  row <- 0L
  for(k in seq_len(n-1L)) {
    children <- hc$merge[k,]
    desc <- integer()
    for(ch in children) {
      child_members <- if(ch<0) -ch else members[[ch]]
      desc <- c(desc,child_members)
      child_height <- if(ch<0) 0 else hc$height[ch]
      row <- row+1L
      lengths[row] <- hc$height[k]-child_height
      incidence[row,child_members] <- TRUE
    }
    members[[k]] <- desc
  }
  lengths[lengths<0 & lengths > -1e-9] <- 0
  if(any(lengths < 0)) stop('Negative dendrogram branch length detected.')
  colnames(incidence) <- species
  list(branches=incidence,lengths=lengths)
}
ptq_tree_beta <- function(branches,lengths,site_species) {
  n <- nrow(site_species)
  occupied <- (branches %*% t(site_species)) > 0
  total <- matrix(0,n,n); turnover <- total; nested <- total
  for(i in seq_len(n)) for(j in seq_len(n)) {
    a <- sum(lengths[occupied[,i] & occupied[,j]])
    b <- sum(lengths[occupied[,i] & !occupied[,j]])
    c <- sum(lengths[!occupied[,i] & occupied[,j]])
    den <- 2*a+b+c
    sor <- if(den>0) (b+c)/den else 0
    sim <- if(a+min(b,c)>0) min(b,c)/(a+min(b,c)) else 0
    total[i,j] <- sor
    turnover[i,j] <- sim
    nested[i,j] <- max(0,sor-sim)
  }
  nm <- rownames(site_species)
  dimnames(total) <- dimnames(turnover) <- dimnames(nested) <- list(nm,nm)
  list(total=total,turnover=turnover,nested=nested)
}
ptq_treebeta_ui <- function() {
  shiny::tabPanel('Trait-tree beta partition',
    shiny::h3('Trait-dendrogram functional beta diversity'),
    shiny::p('Partition pairwise functional dissimilarity into replacement and nestedness-resultant components using a trait dendrogram and presence/absence. This is a tree-based method, not convex-hull functional beta partitioning.'),
    shiny::fluidRow(
      shiny::column(4,shiny::wellPanel(
        shiny::fileInput('tb_traits_file','Species traits CSV',accept='.csv'),
        shiny::fileInput('tb_abundance_file','Site abundance CSV',accept='.csv'),
        shiny::uiOutput('tb_traits_choice'),
        shiny::actionButton('tb_run','Calculate trait-tree beta',class='btn-primary'),
        shiny::downloadButton('tb_example_traits','Download example traits'),
        shiny::downloadButton('tb_example_abundance','Download example abundance'),
        shiny::p('Numeric traits are standardized; average-linkage clustering defines branch lengths. Abundance is converted to presence/absence.'))),
      shiny::column(8,
        shiny::h4('Data validation'),shiny::verbatimTextOutput('tb_validation'),
        shiny::h4('Component to display'),
        shiny::selectInput('tb_component','Beta component',choices=c('Total Sørensen'='total','Turnover (Simpson)'='turnover','Nestedness-resultant'='nested')),
        shiny::h4('Pairwise distance matrix'),shiny::tableOutput('tb_matrix'),
        shiny::downloadButton('tb_matrix_csv','Download selected matrix'),
        shiny::h4('All site-pair components'),shiny::tableOutput('tb_pairs'),
        shiny::downloadButton('tb_pairs_csv','Download all site pairs'),
        shiny::h4('Trait dendrogram'),shiny::plotOutput('tb_tree',height='400px'),
        shiny::h4('Method notes'),shiny::verbatimTextOutput('tb_notes')))
  )
}
ptq_treebeta_server <- function(input,output,session) {
  output$tb_example_traits <- shiny::downloadHandler(
    filename=function()'fd_traits_example.csv',
    content=function(file)file.copy(file.path('data','fd_traits_example.csv'),file,overwrite=TRUE))
  output$tb_example_abundance <- shiny::downloadHandler(
    filename=function()'fd_abundance_example.csv',
    content=function(file)file.copy(file.path('data','fd_abundance_example.csv'),file,overwrite=TRUE))
  getcsv <- function(f) {
    shiny::req(f)
    tryCatch(utils::read.csv(f$datapath,check.names=FALSE,stringsAsFactors=FALSE,
      na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
      error=function(e)shiny::validate(shiny::need(FALSE,paste('CSV error:',conditionMessage(e)))))
  }
  input_data <- shiny::reactive({
    t <- getcsv(input$tb_traits_file); a <- getcsv(input$tb_abundance_file)
    names(t) <- trimws(names(t)); names(a) <- trimws(names(a))
    shiny::validate(shiny::need(!anyDuplicated(names(t)) && !anyDuplicated(names(a)),'Duplicate column names.'),
      shiny::need('species_id' %in% names(t),'Traits require species_id.'),
      shiny::need('site_id' %in% names(a),'Abundance requires site_id.'),
      shiny::need(nrow(a)>=2 && nrow(a)<=250,'Use between 2 and 250 sites for this interactive analysis.'))
    t$species_id <- trimws(as.character(t$species_id))
    a$site_id <- trimws(as.character(a$site_id))
    shiny::validate(shiny::need(!anyNA(t$species_id) && all(nzchar(t$species_id)) &&
      !anyDuplicated(t$species_id),'Invalid or duplicate species_id.'),
      shiny::need(!anyNA(a$site_id) && all(nzchar(a$site_id)) &&
      !anyDuplicated(a$site_id),'Invalid or duplicate site_id.'))
    species <- setdiff(names(a),'site_id')
    shiny::validate(shiny::need(length(species)>=2,'At least two species columns required.'))
    for(nm in species) {
      raw <- a[[nm]]; v <- suppressWarnings(as.numeric(raw))
      shiny::validate(shiny::need(all(is.na(raw) | is.finite(v)),paste('Invalid abundance in',nm)))
      a[[nm]] <- v
    }
    m <- as.matrix(a[species]); storage.mode(m) <- 'double'
    shiny::validate(shiny::need(!anyNA(m) && all(is.finite(m) & m>=0),'Missing, negative or infinite abundance.'),
      shiny::need(all(rowSums(m)>0),'Sites with zero total abundance are not supported.'))
    numeric_traits <- character()
    for(nm in setdiff(names(t),'species_id')) {
      raw <- t[[nm]]; v <- suppressWarnings(as.numeric(raw))
      if(all(is.na(raw) | is.finite(v))) {t[[nm]] <- v; numeric_traits <- c(numeric_traits,nm)}
    }
    shiny::validate(shiny::need(length(numeric_traits)>=1,'At least one numeric trait required.'))
    list(t=t,a=a,m=m,species=species,numeric=numeric_traits)
  })
  output$tb_traits_choice <- shiny::renderUI({
    x <- input_data()
    shiny::selectInput('tb_selected','Select numeric traits',choices=x$numeric,
      selected=x$numeric,multiple=TRUE)
  })
  output$tb_validation <- shiny::renderText({
    if(is.null(input$tb_traits_file) || is.null(input$tb_abundance_file))
      return('Upload the two CSVs or download the included example pair.')
    x <- input_data(); used <- x$species[colSums(x$m)>0]
    missing <- setdiff(used,x$t$species_id)
    paste0('Sites: ',nrow(x$a),' | Observed species: ',length(used),
      '\nSelected traits: ',length(intersect(input$tb_selected,x$numeric)),
      '\nSpecies lacking trait records: ',if(length(missing))paste(missing,collapse=', ') else 'none',
      '\nOnly presence/absence affects this partition.')
  })
  calculated <- shiny::eventReactive(input$tb_run,{
    x <- input_data()
    selected <- intersect(input$tb_selected,x$numeric)
    shiny::validate(shiny::need(length(selected)>=1,'Choose at least one numeric trait.'))
    used <- x$species[colSums(x$m)>0]
    shiny::validate(shiny::need(length(used)>=2,'At least two observed species required.'),
      shiny::need(length(used)<=500,'Limit of 500 observed species for interactive clustering.'),
      shiny::need(all(used %in% x$t$species_id),'Missing trait records for observed species.'))
    t <- x$t[match(used,x$t$species_id),selected,drop=FALSE]
    shiny::validate(shiny::need(all(stats::complete.cases(t)) && all(is.finite(as.matrix(t))),
      'Missing or nonfinite selected traits for observed species.'))
    sdv <- vapply(t,stats::sd,numeric(1))
    shiny::validate(shiny::need(all(is.finite(sdv) & sdv>0),'Remove constant traits.'))
    z <- scale(as.matrix(t))
    rownames(z) <- used
    hc <- stats::hclust(stats::dist(z),method='average')
    branch <- ptq_tree_branch_matrix(hc,used)
    incidence <- (x$m[,match(used,x$species),drop=FALSE]>0)*1
    rownames(incidence) <- x$a$site_id
    beta <- ptq_tree_beta(branch$branches,branch$lengths,incidence)
    pairs <- utils::combn(seq_len(nrow(incidence)),2)
    table <- data.frame(site_1=rownames(incidence)[pairs[1,]],
      site_2=rownames(incidence)[pairs[2,]],
      beta_total=beta$total[t(pairs)],beta_turnover=beta$turnover[t(pairs)],
      beta_nestedness=beta$nested[t(pairs)])
    list(beta=beta,pairs=table,hc=hc,branches=branch)
  })
  selected_matrix <- shiny::reactive({
    z <- calculated()
    z$beta[[input$tb_component]]
  })
  output$tb_matrix <- shiny::renderTable({
    m <- selected_matrix(); as.data.frame(m[seq_len(min(10,nrow(m))),seq_len(min(10,ncol(m))),drop=FALSE])
  },rownames=TRUE,striped=TRUE,digits=5)
  output$tb_pairs <- shiny::renderTable(utils::head(calculated()$pairs,100),striped=TRUE,digits=5)
  output$tb_matrix_csv <- shiny::downloadHandler(
    filename=function()paste0('PlantTraitQuant_trait_tree_',input$tb_component,'.csv'),
    content=function(file)utils::write.csv(selected_matrix(),file,row.names=TRUE))
  output$tb_pairs_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_trait_tree_beta_pairs.csv',
    content=function(file)utils::write.csv(calculated()$pairs,file,row.names=FALSE))
  output$tb_tree <- shiny::renderPlot({
    graphics::par(mar=c(4,4,3,1))
    graphics::plot(calculated()$hc,main='Species trait dendrogram',xlab='',sub='',cex=.8)
  },width=900,height=500,res=96)
  output$tb_notes <- shiny::renderText({
    z <- calculated(); residual <- z$pairs$beta_total-z$pairs$beta_turnover-z$pairs$beta_nestedness
    paste0('Partition: total Sørensen = Simpson turnover + nestedness-resultant component.',
      '\nMaximum numerical decomposition residual: ',signif(max(abs(residual)),4),
      '\nTraits standardized globally; average-linkage Euclidean dendrogram.',
      '\nDendrogram branch lengths weight shared and unique functional history.',
      '\nThis tree-based approach is not equivalent to functional hypervolume overlap.',
      '\nThese are presence/absence components, not abundance-weighted beta diversity.')
  })
}
