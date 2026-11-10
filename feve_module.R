# PlantTraitQuant 2.6.0: Functional evenness from a minimum spanning tree.
# Villéger-style FEve, implemented using base R and Shiny only.
ptq_feve_mst <- function(distance_matrix) {
  n <- nrow(distance_matrix)
  if (n < 2L) stop('Minimum spanning tree needs two or more species.')
  visited <- rep(FALSE,n); visited[1] <- TRUE
  edge_from <- integer(n-1L); edge_to <- integer(n-1L)
  edge_length <- numeric(n-1L)
  for (k in seq_len(n-1L)) {
    candidates <- which(outer(visited,!visited,'&'),arr.ind=TRUE)
    values <- distance_matrix[candidates]
    pick <- which.min(values)
    i <- candidates[pick,1]; j <- candidates[pick,2]
    edge_from[k] <- i; edge_to[k] <- j; edge_length[k] <- values[pick]
    visited[j] <- TRUE
  }
  list(from=edge_from,to=edge_to,length=edge_length)
}
ptq_feve_value <- function(distance_matrix,abundances) {
  idx <- which(abundances>0)
  s <- length(idx)
  if(s<3L) return(NA_real_)
  p <- abundances[idx]/sum(abundances[idx])
  tree <- ptq_feve_mst(distance_matrix[idx,idx,drop=FALSE])
  if (any(!is.finite(tree$length)) || any(tree$length<0)) return(NA_real_)
  ew <- tree$length/(p[tree$from]+p[tree$to])
  if (sum(ew)<=0) return(NA_real_)
  pew <- ew/sum(ew)
  baseline <- 1/(s-1)
  val <- (sum(pmin(pew,baseline))-baseline)/(1-baseline)
  max(0,min(1,val))
}
ptq_feve_ui <- function() {
  shiny::tabPanel('Functional evenness (FEve)',
    shiny::h3('Functional Evenness (FEve)'),
    shiny::p('Site-level FEve from a minimum spanning tree of standardized numeric species traits. Each tab has its own input files.'),
    shiny::fluidRow(
      shiny::column(4,shiny::div(class='well',
        shiny::fileInput('feve_traits_file','1. Species traits CSV',accept='.csv'),
        shiny::fileInput('feve_abund_file','2. Site abundance CSV',accept='.csv'),
        shiny::uiOutput('feve_trait_ui'),
        shiny::actionButton('feve_run','Calculate FEve',class='btn-primary'),
        shiny::downloadButton('feve_traits_example','Download example traits'),
        shiny::downloadButton('feve_abund_example','Download example abundance'),
        shiny::p('Numeric traits are standardized globally. Species distances use Euclidean distance. Each site needs at least three species.'))),
      shiny::column(8,
        shiny::h4('Data validation'),shiny::verbatimTextOutput('feve_validation'),
        shiny::h4('Site-level FEve'),shiny::tableOutput('feve_results'),
        shiny::downloadButton('feve_csv','Download FEve results CSV'),
        shiny::h4('FEve across sites'),shiny::plotOutput('feve_plot',height='400px'),
        shiny::h4('Methods and warnings'),shiny::verbatimTextOutput('feve_notes'))))
}
ptq_feve_server <- function(input,output,session) {
  output$feve_traits_example <- shiny::downloadHandler(
    filename=function()'fd_traits_example.csv',
    content=function(file)file.copy(file.path('data','fd_traits_example.csv'),file,overwrite=TRUE))
  output$feve_abund_example <- shiny::downloadHandler(
    filename=function()'fd_abundance_example.csv',
    content=function(file)file.copy(file.path('data','fd_abundance_example.csv'),file,overwrite=TRUE))
  read_csv <- function(upload) {
    shiny::req(upload)
    tryCatch(utils::read.csv(upload$datapath,check.names=FALSE,
      stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
      error=function(e)shiny::validate(shiny::need(FALSE,paste('CSV read error:',conditionMessage(e)))))
  }
  datasets <- shiny::reactive({
    traits <- read_csv(input$feve_traits_file)
    sites <- read_csv(input$feve_abund_file)
    names(traits) <- trimws(names(traits)); names(sites) <- trimws(names(sites))
    shiny::validate(
      shiny::need(!anyDuplicated(names(traits)) && !anyDuplicated(names(sites)), 'Duplicate column names found.'),
      shiny::need('species_id' %in% names(traits),'Trait file requires species_id.'),
      shiny::need('site_id' %in% names(sites),'Abundance file requires site_id.'),
      shiny::need(nrow(traits)>=3 && nrow(sites)>0,'At least three species and one site required.'))
    traits$species_id <- trimws(as.character(traits$species_id))
    sites$site_id <- trimws(as.character(sites$site_id))
    shiny::validate(
      shiny::need(all(!is.na(traits$species_id) & nzchar(traits$species_id)) &&
        !anyDuplicated(traits$species_id),'Species identifiers must be nonempty and unique.'),
      shiny::need(all(!is.na(sites$site_id) & nzchar(sites$site_id)) &&
        !anyDuplicated(sites$site_id),'Site identifiers must be nonempty and unique.'))
    species_cols <- setdiff(names(sites),'site_id')
    trait_cols <- setdiff(names(traits),'species_id')
    shiny::validate(shiny::need(length(species_cols)>0,'No abundance species columns found.'),
      shiny::need(length(trait_cols)>=2,'At least two numeric trait columns are needed.'))
    for(nm in species_cols) {
      raw <- sites[[nm]]; num <- suppressWarnings(as.numeric(raw))
      shiny::validate(shiny::need(all(is.na(raw) | is.finite(num)),paste('Invalid abundance in',nm)))
      sites[[nm]] <- num
    }
    m <- as.matrix(sites[species_cols])
    shiny::validate(shiny::need(!anyNA(m),'Missing abundance cells: please correct the CSV.'),
      shiny::need(all(is.finite(m) & m>=0),'Abundances must be finite and nonnegative.'))
    for(nm in trait_cols) {
      if(is.character(traits[[nm]])) {
        raw <- traits[[nm]]; num <- suppressWarnings(as.numeric(raw))
        if(all(is.na(raw) | is.finite(num))) traits[[nm]] <- num
      }
    }
    numeric_cols <- trait_cols[vapply(traits[trait_cols],is.numeric,logical(1))]
    shiny::validate(shiny::need(length(numeric_cols)>=2,'At least two numeric traits required.'))
    list(traits=traits,sites=sites,m=m,species_cols=species_cols,numeric_cols=numeric_cols)
  })
  output$feve_trait_ui <- shiny::renderUI({
    d <- datasets()
    shiny::selectInput('feve_selected_traits','3. Select numerical traits',
      choices=d$numeric_cols,selected=d$numeric_cols,multiple=TRUE)
  })
  output$feve_validation <- shiny::renderText({
    if(is.null(input$feve_traits_file) || is.null(input$feve_abund_file))
      return('Upload both CSV files or use the included examples.')
    d <- datasets(); selected <- intersect(input$feve_selected_traits,d$numeric_cols)
    active <- d$species_cols[colSums(d$m)>0]
    unmatched <- setdiff(active,d$traits$species_id)
    paste0('Sites: ',nrow(d$sites),' | Species trait records: ',nrow(d$traits),
      '\nSelected numerical traits: ',length(selected),
      '\nPresent species missing from traits: ',
      if(length(unmatched))paste(unmatched,collapse=', ') else 'none',
      '\nMissing trait cells (selected columns): ',
      if(length(selected))sum(is.na(d$traits[selected])) else 0,
      '\nSites with fewer than three present species: ',sum(rowSums(d$m>0)<3))
  })
  computed <- shiny::eventReactive(input$feve_run,{
    d <- datasets(); selected <- intersect(input$feve_selected_traits,d$numeric_cols)
    shiny::validate(shiny::need(length(selected)>=2,'Select at least two numeric traits.'))
    active <- d$species_cols[colSums(d$m)>0]
    shiny::validate(shiny::need(length(active)>=3,'At least three species must be present across sites.'),
      shiny::need(length(active)<=1200,'For this interactive FEve module, limit to 1200 species per run.'))
    unmatched <- setdiff(active,d$traits$species_id)
    shiny::validate(shiny::need(!length(unmatched),
      paste('No trait records for:',paste(unmatched,collapse=', '))))
    t <- d$traits[match(active,d$traits$species_id),selected,drop=FALSE]
    shiny::validate(shiny::need(!anyNA(t) && all(is.finite(as.matrix(t))),
      'Selected traits must be complete and finite for species occurring in the community.'))
    sds <- vapply(t,stats::sd,numeric(1))
    shiny::validate(shiny::need(all(is.finite(sds) & sds>0),
      'At least one selected trait is constant; remove constant traits.'))
    distance <- as.matrix(stats::dist(scale(as.matrix(t)),method='euclidean'))
    m <- d$m[,match(active,d$species_cols),drop=FALSE]
    values <- vapply(seq_len(nrow(m)),function(i)ptq_feve_value(distance,m[i,]),numeric(1))
    out <- data.frame(site_id=d$sites$site_id,
      species_richness=rowSums(m>0),FEve=values,check.names=FALSE)
    list(table=out)
  })
  output$feve_results <- shiny::renderTable(computed()$table,striped=TRUE,digits=5)
  output$feve_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_FEve.csv',
    content=function(file)utils::write.csv(computed()$table,file,row.names=FALSE,na=''))
  output$feve_plot <- shiny::renderPlot({
    d <- computed()$table
    graphics::par(mar=c(7,4.5,2,1))
    graphics::barplot(d$FEve,names.arg=d$site_id,las=2,col='#19745b',
      ylab='Functional Evenness (FEve)',ylim=c(0,1),main='FEve by site')
  },width=850,height=510,res=96)
  output$feve_notes <- shiny::renderText({
    d <- computed()$table
    paste0('Minimum spanning tree (MST) of standardized Euclidean trait distances.',
      '\nEdge weight EW = edge distance / (relative abundance of its two species).',
      '\nNormalized edge weights are capped at 1/(S-1) to calculate FEve.',
      '\nFEve ranges from 0 to 1; it is undefined for sites with fewer than 3 species',
      ' or zero total MST edge length.',
      '\nUndefined sites: ',sum(is.na(d$FEve)),
      '\nResults should be independently compared against published reference implementations.')
  })
}
