
# Multidimensional convex-hull FRic (not a probabilistic hypervolume).
ptq_hull_ui <- function() {
 shiny::tabPanel('Multidimensional FRic',
  shiny::h3('Multidimensional functional richness and trait-space volume'),
  shiny::p('Calculate convex-hull volume in 2 to 5 standardized numerical trait dimensions using geometry::convhulln. Upload your own trait and abundance files.'),
  shiny::fluidRow(
   shiny::column(4,shiny::div(class='well',
    shiny::fileInput('hv_traits_file','Species traits CSV',accept='.csv'),
    shiny::fileInput('hv_abundance_file','Site abundance CSV',accept='.csv'),
    shiny::uiOutput('hv_choices'),
    shiny::numericInput('hv_dim','Number of dimensions',2,min=2,max=5,step=1),
    shiny::actionButton('hv_run','Calculate multidimensional FRic',class='btn-primary'),
    shiny::downloadButton('hv_traits_example','Download example traits'),
    shiny::downloadButton('hv_abundance_example','Download example abundance'))),
   shiny::column(8,
    shiny::h4('Data validation'),shiny::verbatimTextOutput('hv_checks'),
    shiny::h4('Site-level convex-hull functional richness'),shiny::tableOutput('hv_results'),
    shiny::downloadButton('hv_results_csv','Download FRic results'),
    shiny::h4('Reference test: unit simplex'),shiny::verbatimTextOutput('hv_test'),
    shiny::h4('Scope and interpretation'),shiny::verbatimTextOutput('hv_notes'))))
}
ptq_hull_server <- function(input,output,session) {
 for(nm in c('traits','abundance'))local({
  key <- nm
  output[[paste0('hv_',key,'_example')]] <- shiny::downloadHandler(
   filename=function()paste0('functional_structure_',key,'.csv'),
   content=function(file)file.copy(file.path('data',paste0('functional_structure_',key,'.csv')),file,overwrite=TRUE))
 })
 read_file <- function(upload) {
  shiny::req(upload)
  x <- tryCatch(utils::read.csv(upload$datapath,check.names=FALSE,
   stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
   error=function(e)shiny::validate(shiny::need(FALSE,conditionMessage(e))))
  names(x) <- trimws(names(x))
  shiny::validate(shiny::need(!anyDuplicated(names(x)),'Duplicate column names.'))
  x
 }
 dat <- shiny::reactive({
  t <- read_file(input$hv_traits_file);a <- read_file(input$hv_abundance_file)
  shiny::validate(shiny::need('species_id'%in%names(t),'Trait CSV needs species_id.'),
   shiny::need('site_id'%in%names(a),'Abundance CSV needs site_id.'))
  t$species_id <- trimws(as.character(t$species_id))
  a$site_id <- trimws(as.character(a$site_id))
  shiny::validate(shiny::need(!anyNA(t$species_id)&&all(nzchar(t$species_id))&&
   !anyDuplicated(t$species_id),'Species IDs must be unique and nonempty.'),
   shiny::need(!anyNA(a$site_id)&&all(nzchar(a$site_id))&&
   !anyDuplicated(a$site_id),'Site IDs must be unique and nonempty.'))
  species <- setdiff(names(a),'site_id')
  shiny::validate(shiny::need(length(species)>=3,'At least three species columns required.'),
   shiny::need(all(species%in%t$species_id),'Some abundance species are absent from traits.'))
  abundance <- sapply(a[species],function(v)suppressWarnings(as.numeric(as.character(v))))
  abundance <- matrix(abundance,nrow=nrow(a),dimnames=list(a$site_id,species))
  shiny::validate(shiny::need(all(is.finite(abundance))&&all(abundance>=0),
   'Abundance must contain finite nonnegative values.'))
  active <- colSums(abundance)>0
  abundance <- abundance[,active,drop=FALSE]
  shiny::validate(shiny::need(ncol(abundance)>=3,'At least three species must occur.'))
  t <- t[match(colnames(abundance),t$species_id),,drop=FALSE]
  possible <- setdiff(names(t),'species_id')
  possible <- possible[vapply(t[possible],function(v){
   z <- suppressWarnings(as.numeric(as.character(v)))
   all(!is.na(z) & is.finite(z))
  },logical(1))]
  list(t=t,a=abundance,possible=possible)
 })
 output$hv_choices <- shiny::renderUI({
  x <- dat()
  shiny::selectInput('hv_traits','Numeric traits',choices=x$possible,
   selected=x$possible,multiple=TRUE)
 })
 output$hv_checks <- shiny::renderText({
  if(is.null(input$hv_traits_file)||is.null(input$hv_abundance_file))
   return('Upload both CSV files or use the examples.')
  x <- dat()
  paste0('Sites: ',nrow(x$a),'\nOccurring species: ',ncol(x$a),
   '\nComplete numerical traits: ',length(x$possible),
   '\ngeometry installed: ',requireNamespace('geometry',quietly=TRUE))
 })
 hull_volume <- function(coords) {
  # Requires at least d+1 points in full affine rank.
  d <- ncol(coords)
  if(nrow(coords)<d+1)return(NA_real_)
  if(qr(sweep(coords,2,colMeans(coords),'-'))$rank<d)return(0)
  tryCatch(as.numeric(geometry::convhulln(coords,options='FA',output.options=TRUE)$vol),
   error=function(e)NA_real_)
 }
 result <- shiny::eventReactive(input$hv_run,{
  shiny::validate(shiny::need(requireNamespace('geometry',quietly=TRUE),
   'Install the geometry package: install.packages("geometry").'))
  x <- dat();sel <- intersect(input$hv_traits,x$possible)
  k <- as.integer(input$hv_dim)
  shiny::validate(shiny::need(is.finite(k)&&k>=2&&k<=5,'Choose 2 to 5 dimensions.'),
   shiny::need(length(sel)>=k,'Select at least as many traits as dimensions.'),
   shiny::need(ncol(x$a)>=k+1,'Need at least dimension+1 occurring species.'))
  traits <- sapply(x$t[sel],function(v)as.numeric(as.character(v)))
  traits <- matrix(traits,nrow=ncol(x$a),dimnames=list(colnames(x$a),sel))
  shiny::validate(shiny::need(all(apply(traits,2,stats::sd)>0),
   'Remove constant traits.'))
  standardized <- scale(traits)
  pc <- stats::prcomp(standardized,center=FALSE,scale.=FALSE)
  shiny::validate(shiny::need(sum(pc$sdev>1e-9)>=k,
   'The trait space has insufficient independent dimensions.'))
  coords <- pc$x[,seq_len(k),drop=FALSE]
  rows <- lapply(seq_len(nrow(x$a)),function(i){
   present <- which(x$a[i,]>0)
   data.frame(site_id=rownames(x$a)[i],species_richness=length(present),
    dimensions=k,
    hull_volume=if(length(present)<k+1)NA_real_ else hull_volume(coords[present,,drop=FALSE]),
    status=if(length(present)<k+1)'Insufficient species' else
     if(qr(sweep(coords[present,,drop=FALSE],2,colMeans(coords[present,,drop=FALSE]),'-'))$rank<k)
      'Degenerate hull (zero volume)' else 'Calculated')
  })
  list(table=do.call(rbind,rows),k=k,variance=100*sum(pc$sdev[seq_len(k)]^2)/sum(pc$sdev^2))
 })
 output$hv_results <- shiny::renderTable(result()$table,striped=TRUE,digits=5)
 output$hv_test <- shiny::renderText({
  if(!requireNamespace('geometry',quietly=TRUE))return('Install geometry to run reference test.')
  # A d-dimensional unit simplex has volume 1/d!.
  d <- 3L
  simplex <- rbind(rep(0,d),diag(d))
  actual <- hull_volume(simplex);expected <- 1/factorial(d)
  paste0('3D unit simplex: expected volume ',format(expected,digits=8),
   ', calculated ',format(actual,digits=8),
   '\n',if(is.finite(actual)&&abs(actual-expected)<1e-7)'PASS' else 'FAIL')
 })
 output$hv_notes <- shiny::renderText({
  z <- result()
  paste0('Dimensions: ',z$k,'; variance represented by chosen PCA axes: ',
   round(z$variance,2),'%.\n',
   'FRic is the convex-hull volume of occurring species in the first selected PCA axes of globally standardized numeric traits.\n',
   'Species abundances determine presence only; hull volume is abundance-insensitive.\n',
   'A site with fewer than dimensions+1 species is NA; affinely degenerate hulls have zero volume.\n',
   'Results depend on selected traits, PCA dimensions and the global species pool.\n',
   'This is multidimensional convex-hull FRic, NOT a probabilistic kernel hypervolume or TPD.\n',
   'The unit-simplex reference test verifies geometric volume computation, not ecological index equivalence.')
 })
 output$hv_results_csv <- shiny::downloadHandler(
  filename=function()'PlantTraitQuant_multidimensional_FRic.csv',
  content=function(file)utils::write.csv(result()$table,file,row.names=FALSE,na=''))
}
