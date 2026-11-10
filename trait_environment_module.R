
# Three-table trait-environment analysis: weighted cross-covariance RLQ-style
# ordination and fourth-corner weighted correlations, with site-label permutations.
ptq_te_ui <- function() {
 shiny::tabPanel('Trait-environment associations',
  shiny::h3('Trait-environment associations and RLQ-style ordination'),
  shiny::p('Three linked CSVs: L = site-by-species abundance, R = site environmental variables, Q = species traits. Numerical variables only in this first implementation.'),
  shiny::fluidRow(
   shiny::column(4,shiny::div(class='well',
    shiny::fileInput('te_l','L: Site abundance CSV',accept='.csv'),
    shiny::fileInput('te_r','R: Site environment CSV',accept='.csv'),
    shiny::fileInput('te_q','Q: Species traits CSV',accept='.csv'),
    shiny::uiOutput('te_vars'),
    shiny::numericInput('te_perm','Site-label permutations',199,min=19,max=1999,step=1),
    shiny::numericInput('te_seed','Random seed',123,min=1,max=1000000),
    shiny::actionButton('te_go','Calculate associations',class='btn-primary'),
    shiny::downloadButton('te_l_example','Example abundance'),
    shiny::downloadButton('te_r_example','Example environment'),
    shiny::downloadButton('te_q_example','Example traits'))),
   shiny::column(8,
    shiny::h4('Data alignment'),shiny::verbatimTextOutput('te_checks'),
    shiny::h4('Weighted fourth-corner correlations'),shiny::tableOutput('te_cor'),
    shiny::downloadButton('te_cor_csv','Download fourth-corner results'),
    shiny::h4('RLQ-style weighted cross-covariance ordination'),shiny::tableOutput('te_axes'),
    shiny::plotOutput('te_plot',height='420px'),
    shiny::downloadButton('te_axes_csv','Download ordination loadings'),
    shiny::h4('Interpretation and scope'),shiny::verbatimTextOutput('te_notes'))))
}
ptq_te_server <- function(input,output,session) {
 for(nm in c('l','r','q')) {
  local({
   name <- nm
   output[[paste0('te_',name,'_example')]] <- shiny::downloadHandler(
    filename=function()paste0('trait_environment_',name,'.csv'),
    content=function(file)file.copy(file.path('data',paste0('trait_environment_',name,'.csv')),file,overwrite=TRUE))
  })
 }
 read_table <- function(upload,id) {
  shiny::req(upload)
  x <- tryCatch(utils::read.csv(upload$datapath,check.names=FALSE,
   stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
   error=function(e)shiny::validate(shiny::need(FALSE,conditionMessage(e))))
  names(x) <- trimws(names(x))
  shiny::validate(shiny::need(id%in%names(x),paste('Missing',id)),
   shiny::need(!anyDuplicated(names(x)),'Duplicate column names.'))
  x[[id]] <- trimws(as.character(x[[id]]))
  shiny::validate(shiny::need(!anyNA(x[[id]]) && all(nzchar(x[[id]])) &&
   !anyDuplicated(x[[id]]),paste('Invalid or duplicate',id)))
  x
 }
 input_tables <- shiny::reactive({
  l <- read_table(input$te_l,'site_id')
  r <- read_table(input$te_r,'site_id')
  q <- read_table(input$te_q,'species_id')
  species <- setdiff(names(l),'site_id')
  shiny::validate(shiny::need(length(species)>=2,'L needs at least two species columns.'),
   shiny::need(all(species%in%q$species_id),
    paste('Missing trait species:',paste(setdiff(species,q$species_id),collapse=', '))),
   shiny::need(all(l$site_id%in%r$site_id),
    paste('Missing environment sites:',paste(setdiff(l$site_id,r$site_id),collapse=', '))),
   shiny::need(nrow(l)>=4,'At least four sites are required.'))
  a <- as.matrix(l[,species,drop=FALSE])
  suppressWarnings(storage.mode(a) <- 'double')
  shiny::validate(shiny::need(all(is.finite(a)) && all(a>=0),
   'L must contain nonnegative finite numeric abundances without missing cells.'),
   shiny::need(sum(a)>0,'Abundance matrix must contain positive values.'))
  r <- r[match(l$site_id,r$site_id),,drop=FALSE]
  q <- q[match(species,q$species_id),,drop=FALSE]
  list(L=a,R=r,Q=q,sites=l$site_id,species=species)
 })
 output$te_vars <- shiny::renderUI({
  x <- input_tables()
  numeric_cols <- function(d,id) {
   names0 <- setdiff(names(d),id)
   names0[vapply(d[names0],function(v){
    z <- suppressWarnings(as.numeric(as.character(v)))
    all(is.na(v)|is.finite(z))
   },logical(1))]
  }
  shiny::tagList(
   shiny::selectInput('te_env','Numeric environmental variables',
    choices=numeric_cols(x$R,'site_id'),
    selected=numeric_cols(x$R,'site_id'),multiple=TRUE),
   shiny::selectInput('te_trait','Numeric species traits',
    choices=numeric_cols(x$Q,'species_id'),
    selected=numeric_cols(x$Q,'species_id'),multiple=TRUE))
 })
 output$te_checks <- shiny::renderText({
  if(is.null(input$te_l)||is.null(input$te_r)||is.null(input$te_q))
   return('Upload all three CSVs or use the bundled example files.')
  x <- input_tables()
  paste0('Aligned sites: ',nrow(x$L),'\nAligned species: ',ncol(x$L),
   '\nSpecies absent from L but present in Q: ',
   length(setdiff(x$Q$species_id,x$species)),
   '\nSelected environment variables: ',length(input$te_env),
   '\nSelected traits: ',length(input$te_trait))
 })
 result <- shiny::eventReactive(input$te_go,{
  x <- input_tables()
  ev <- intersect(input$te_env,setdiff(names(x$R),'site_id'))
  tr <- intersect(input$te_trait,setdiff(names(x$Q),'species_id'))
  shiny::validate(shiny::need(length(ev)>0 && length(tr)>0,
   'Select at least one environmental variable and one trait.'))
  R <- sapply(x$R[ev],function(z)as.numeric(as.character(z)))
  Q <- sapply(x$Q[tr],function(z)as.numeric(as.character(z)))
  R <- matrix(R,nrow=nrow(x$L),dimnames=list(x$sites,ev))
  Q <- matrix(Q,nrow=ncol(x$L),dimnames=list(x$species,tr))
  shiny::validate(shiny::need(all(is.finite(R)) && all(is.finite(Q)),
   'Selected R and Q columns must be complete finite numeric values.'))
  L <- x$L
  active_sites <- rowSums(L)>0
  active_species <- colSums(L)>0
  shiny::validate(shiny::need(sum(active_sites)>=4 && sum(active_species)>=2,
   'Need at least four nonempty sites and two occurring species.'))
  L <- L[active_sites,active_species,drop=FALSE]
  R <- R[active_sites,,drop=FALSE]
  Q <- Q[active_species,,drop=FALSE]
  wr <- rowSums(L)/sum(L);wc <- colSums(L)/sum(L)
  standardize <- function(a,w) {
   center <- colSums(sweep(a,1,w,'*'))
   b <- sweep(a,2,center,'-')
   sd <- sqrt(colSums(sweep(b*b,1,w,'*')))
   shiny::validate(shiny::need(all(is.finite(sd) & sd>1e-12),
    paste('Constant or invalid variable:',paste(colnames(a)[!is.finite(sd)|sd<=1e-12],collapse=', '))))
   sweep(b,2,sd,'/')
  }
  E <- standardize(R,wr)
  T <- standardize(Q,wc)
  P <- L/sum(L)
  # Weighted fourth-corner association C = E' P T
  calc_cross <- function(e) t(e)%*%P%*%T
  cross <- calc_cross(E)
  perms <- as.integer(input$te_perm);seed <- as.integer(input$te_seed)
  shiny::validate(shiny::need(is.finite(perms) && perms>=19 && perms<=1999,
   'Permutations must be between 19 and 1999.'),
   shiny::need(is.finite(seed) && seed>=1,'Seed must be positive.'))
  set.seed(seed)
  exceed <- matrix(0L,nrow=nrow(cross),ncol=ncol(cross))
  for(b in seq_len(perms)) {
   perm <- calc_cross(E[sample.int(nrow(E)),,drop=FALSE])
   exceed <- exceed+(abs(perm)>=abs(cross)-1e-12)
  }
  p <- (exceed+1)/(perms+1)
  # BH correction across all environment-trait pairs
  padj <- matrix(stats::p.adjust(as.vector(p),method='BH'),nrow=nrow(p))
  tab <- expand.grid(environment=rownames(cross),trait=colnames(cross),
   stringsAsFactors=FALSE)
  tab$weighted_association <- as.vector(cross)
  tab$site_permutation_p <- as.vector(p)
  tab$BH_adjusted_p <- as.vector(padj)
  # Singular vectors of weighted standardized cross-covariance
  sv <- svd(cross,nu=min(dim(cross)),nv=min(dim(cross)))
  naxes <- length(sv$d)
  axes <- rbind(
   data.frame(type='Environment',variable=rownames(cross),
    axis1=sv$u[,1],axis2=if(naxes>=2)sv$u[,2] else NA_real_),
   data.frame(type='Trait',variable=colnames(cross),
    axis1=sv$v[,1],axis2=if(naxes>=2)sv$v[,2] else NA_real_))
  list(tab=tab,axes=axes,values=sv$d,nsite=nrow(L),nspecies=ncol(L),perms=perms)
 })
 output$te_cor <- shiny::renderTable(result()$tab,striped=TRUE,digits=4)
 output$te_axes <- shiny::renderTable(result()$axes,striped=TRUE,digits=4)
 output$te_plot <- shiny::renderPlot({
  z <- result()
  a <- z$axes
  if(all(is.na(a$axis2))) {
   graphics::plot.new();graphics::text(0.5,0.5,'Two axes require at least two traits and two environment variables.')
  } else {
   graphics::plot(a$axis1,a$axis2,type='n',xlab='Cross-covariance axis 1',
    ylab='Cross-covariance axis 2',main='RLQ-style trait-environment loadings')
   graphics::abline(h=0,v=0,col='grey80')
   graphics::points(a$axis1,a$axis2,pch=ifelse(a$type=='Environment',16,17),
    col=ifelse(a$type=='Environment','steelblue4','darkgreen'))
   graphics::text(a$axis1,a$axis2,labels=a$variable,pos=3,cex=.8)
   graphics::legend('topright',legend=c('Environment','Trait'),
    pch=c(16,17),col=c('steelblue4','darkgreen'),bty='n')
  }
 })
 output$te_notes <- shiny::renderText({
  z <- result()
  paste0('Active sites: ',z$nsite,'; occurring species: ',z$nspecies,
   '; site permutations: ',z$perms,'.\n',
   'Associations are abundance-weighted cross-moments of standardized R and Q, not ordinary Pearson correlations.\n',
   'Permutations shuffle site-environment rows while keeping abundance and traits fixed; this tests a site-label null only.\n',
   'BH-adjusted p-values correct across all displayed environment-trait pairs.\n',
   'Ordination is an SVD of the weighted cross-covariance matrix (RLQ-style), NOT a full ade4::rlq implementation.\n',
   'No model-2/4 combined permutation test, spatial blocking or covariate adjustment is implemented.\n',
   'For confirmatory fourth-corner inference use appropriate model-2/model-4 testing in ade4, especially with spatial dependence.')
 })
 output$te_cor_csv <- shiny::downloadHandler(
  filename=function()'PlantTraitQuant_fourth_corner_associations.csv',
  content=function(file)utils::write.csv(result()$tab,file,row.names=FALSE))
 output$te_axes_csv <- shiny::downloadHandler(
  filename=function()'PlantTraitQuant_RLQ_style_loadings.csv',
  content=function(file)utils::write.csv(result()$axes,file,row.names=FALSE))
}
