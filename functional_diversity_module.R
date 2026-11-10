# Community functional diversity: standardized Euclidean FDis and RaoQ.
# Requires only shiny and base R.
ptq_fd_ui <- function() {
  shiny::tabPanel('Community functional diversity',
    shiny::h3('Community functional diversity: FDis and RaoQ'),
    shiny::p('Upload a species traits CSV and a site-by-species abundance CSV. Species names must match exactly. This first release supports numeric traits.'),
    shiny::fluidRow(
      shiny::column(4, shiny::div(class='well',
        shiny::fileInput('fd_traits_file','1. Upload species traits CSV',accept='.csv'),
        shiny::fileInput('fd_abundance_file','2. Upload site abundance CSV',accept='.csv'),
        shiny::uiOutput('fd_trait_selection'),
        shiny::selectInput('fd_weight','Abundance weighting',choices=c('Use relative abundance'='abundance','Presence/absence (equal species weight)'='presence')),
        shiny::actionButton('fd_run','Calculate FDis and RaoQ',class='btn-primary'),
        shiny::downloadButton('fd_example_traits','Example traits CSV'),
        shiny::downloadButton('fd_example_abundance','Example abundance CSV'),
        shiny::p('Missing trait observations, unmatched species with positive abundance, and missing abundance values stop the calculation.'))),
      shiny::column(8,
        shiny::h4('Input checks'),shiny::verbatimTextOutput('fd_checks'),
        shiny::h4('Site-level functional diversity'),shiny::tableOutput('fd_results'),
        shiny::downloadButton('fd_results_csv','Download functional diversity CSV'),
        shiny::h4('Site comparison'),shiny::plotOutput('fd_plot',height='390px'),
        shiny::p('FDis is the abundance-weighted mean Euclidean distance to the weighted centroid in standardized trait space. RaoQ is the sum of p_i p_j times standardized Euclidean species distance. Values depend on selected traits and their scaling.'))))
}
ptq_fd_server <- function(input,output,session) {
  output$fd_example_traits <- shiny::downloadHandler(
    filename=function() 'fd_traits_example.csv',
    content=function(file) file.copy(file.path('data','fd_traits_example.csv'),file,overwrite=TRUE))
  output$fd_example_abundance <- shiny::downloadHandler(
    filename=function() 'fd_abundance_example.csv',
    content=function(file) file.copy(file.path('data','fd_abundance_example.csv'),file,overwrite=TRUE))
  read_csv <- function(file) {
    shiny::req(file)
    tryCatch(utils::read.csv(file$datapath,check.names=FALSE,stringsAsFactors=FALSE,
                             na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
      error=function(e) shiny::validate(shiny::need(FALSE,paste('Cannot read CSV:',conditionMessage(e)))))
  }
  files <- shiny::reactive({
    t <- read_csv(input$fd_traits_file)
    a <- read_csv(input$fd_abundance_file)
    names(t) <- trimws(names(t)); names(a) <- trimws(names(a))
    shiny::validate(
      shiny::need(!anyDuplicated(names(t)) && !anyDuplicated(names(a)),'Duplicate column names detected.'),
      shiny::need('species_id' %in% names(t),'Traits CSV needs species_id.'),
      shiny::need('site_id' %in% names(a),'Abundance CSV needs site_id.'),
      shiny::need(nrow(t)>=2 && nrow(a)>=1,'At least two trait species and one site are required.'))
    t$species_id <- trimws(as.character(t$species_id))
    a$site_id <- trimws(as.character(a$site_id))
    shiny::validate(
      shiny::need(!anyNA(t$species_id) && all(nzchar(t$species_id)) && !anyDuplicated(t$species_id),'species_id must be unique and nonblank.'),
      shiny::need(!anyNA(a$site_id) && all(nzchar(a$site_id)) && !anyDuplicated(a$site_id),'site_id must be unique and nonblank.'))
    tc <- setdiff(names(t),'species_id'); ac <- setdiff(names(a),'site_id')
    shiny::validate(shiny::need(length(tc)>0 && length(ac)>0,'Both files need data columns.'))
    for(nm in tc) {
      v <- suppressWarnings(as.numeric(t[[nm]]))
      if(all(is.na(t[[nm]]) | is.finite(v))) t[[nm]] <- v
    }
    for(nm in ac) {
      v <- suppressWarnings(as.numeric(a[[nm]]))
      shiny::validate(shiny::need(all(is.na(a[[nm]]) | is.finite(v)),paste('Non-numeric abundance:',nm)))
      a[[nm]] <- v
    }
    traits <- tc[vapply(t[tc],is.numeric,logical(1))]
    shiny::validate(shiny::need(length(traits)>=1,'At least one numeric trait required.'),
      shiny::need(!any(as.matrix(a[ac])<0,na.rm=TRUE),'Negative abundances are invalid.'),
      shiny::need(!anyNA(a[ac]),'Missing abundance values found. Correct CSV before calculating.'),
      shiny::need(all(is.finite(as.matrix(a[ac]))),'Nonfinite abundance values found.'))
    list(traits=t,abundance=a,trait_names=traits,abundance_names=ac)
  })
  output$fd_trait_selection <- shiny::renderUI({
    x <- files()
    shiny::selectInput('fd_traits','3. Select numeric traits',choices=x$trait_names,
                       selected=x$trait_names,multiple=TRUE)
  })
  output$fd_checks <- shiny::renderText({
    if(is.null(input$fd_traits_file) || is.null(input$fd_abundance_file))
      return('Upload both CSV files or download the paired examples.')
    x <- files()
    present <- x$abundance_names[colSums(x$abundance[x$abundance_names])>0]
    missing <- setdiff(present,x$traits$species_id)
    selected <- intersect(input$fd_traits,x$trait_names)
    paste0('Sites: ',nrow(x$abundance),' | Trait species: ',nrow(x$traits),
      '\nSpecies columns: ',length(x$abundance_names),
      '\nSelected numeric traits: ',length(selected),
      '\nSpecies with positive abundance missing trait records: ',
      if(length(missing))paste(missing,collapse=', ') else 'none',
      '\nMissing selected trait cells: ',if(length(selected))sum(is.na(x$traits[selected])) else 0,
      '\nZero-total sites: ',sum(rowSums(x$abundance[x$abundance_names])==0),
      '\nTrait scaling: mean 0, standard deviation 1 across observed species.')
  })
  calculated <- shiny::eventReactive(input$fd_run,{
    x <- files()
    selected <- intersect(input$fd_traits,x$trait_names)
    shiny::validate(shiny::need(length(selected)>=1,'Select at least one numeric trait.'))
    abundance <- as.matrix(x$abundance[x$abundance_names])
    present <- x$abundance_names[colSums(abundance)>0]
    missing <- setdiff(present,x$traits$species_id)
    shiny::validate(shiny::need(!length(missing),paste('Positive-abundance species without trait records:',paste(missing,collapse=', '))))
    used <- x$traits[match(present,x$traits$species_id),selected,drop=FALSE]
    shiny::validate(shiny::need(!anyNA(used),'Missing numeric trait values for species present in sites.'),
      shiny::need(all(is.finite(as.matrix(used))),'Nonfinite numeric traits found.'),
      shiny::need(length(present)>=2,'At least two species with positive abundance are required.'))
    standard <- apply(used,2,stats::sd)
    shiny::validate(shiny::need(all(is.finite(standard) & standard>0),
      'Selected traits must vary across species with positive abundance. Remove constant traits.'))
    z <- scale(as.matrix(used),center=TRUE,scale=TRUE)
    rownames(z) <- present
    dm <- as.matrix(stats::dist(z,method='euclidean'))
    res <- lapply(seq_len(nrow(abundance)),function(i){
      raw <- abundance[i,match(present,x$abundance_names)]
      if(identical(input$fd_weight,'presence')) raw <- as.numeric(raw>0)
      n <- sum(raw)
      if(n==0) return(c(species_richness=0,FDis=NA_real_,RaoQ=NA_real_))
      p <- raw/n
      centroid <- colSums(z*p)
      radial <- sqrt(rowSums(sweep(z,2,centroid,'-')^2))
      c(species_richness=sum(raw>0),FDis=sum(p*radial),
        RaoQ=as.numeric(crossprod(p,dm%*%p)))
    })
    results <- data.frame(site_id=x$abundance$site_id,do.call(rbind,res),row.names=NULL)
    results
  })
  output$fd_results <- shiny::renderTable(calculated(),striped=TRUE,digits=4)
  output$fd_results_csv <- shiny::downloadHandler(
    filename=function() 'PlantTraitQuant_FDis_RaoQ.csv',
    content=function(file) utils::write.csv(calculated(),file,row.names=FALSE,na=''))
  output$fd_plot <- shiny::renderPlot({
    d <- calculated()
    if(!any(is.finite(d$FDis))) return(invisible(NULL))
    op <- graphics::par(mar=c(7,4,2,1))
    on.exit(graphics::par(op))
    graphics::barplot(d$FDis,names.arg=d$site_id,las=2,col='#218765',
                      ylab='Functional dispersion (FDis)',main='Site-level functional dispersion')
  },width=850,height=420,res=96)
}
