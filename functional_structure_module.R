
# Exploratory functional structure: redundancy, originality, specialization,
# vulnerability, null-model SES, pooled scale curves, and 2D hull proxy.
ptq_structure_ui <- function() {
 shiny::tabPanel('Functional structure & scales',
  shiny::h3('Functional structure, null models and spatial-scale summaries'),
  shiny::p('Upload your own paired species traits and site abundance CSVs. Optional site grouping CSV supports pooled comparisons across defined spatial units.'),
  shiny::fluidRow(
   shiny::column(4,shiny::div(class='well',
    shiny::fileInput('fs_traits_file','Species trait CSV',accept='.csv'),
    shiny::fileInput('fs_abundance_file','Site abundance CSV',accept='.csv'),
    shiny::fileInput('fs_groups_file','Optional site groups CSV',accept='.csv'),
    shiny::uiOutput('fs_trait_choices'),
    shiny::numericInput('fs_reps','Null permutations',199,min=19,max=999,step=1),
    shiny::numericInput('fs_seed','Random seed',123,min=1,max=1000000),
    shiny::actionButton('fs_run','Calculate functional structure',class='btn-primary'),
    shiny::downloadButton('fs_traits_example','Example traits'),
    shiny::downloadButton('fs_abundance_example','Example abundance'),
    shiny::downloadButton('fs_groups_example','Example site groups'))),
   shiny::column(8,
    shiny::h4('Validation'),shiny::verbatimTextOutput('fs_checks'),
    shiny::h4('Site functional structure and null-model SES'),shiny::tableOutput('fs_site'),
    shiny::downloadButton('fs_site_csv','Download site indices'),
    shiny::h4('Species originality and specialization'),shiny::tableOutput('fs_species'),
    shiny::downloadButton('fs_species_csv','Download species indices'),
    shiny::h4('Pooled site accumulation and group summaries'),shiny::tableOutput('fs_scale'),
    shiny::downloadButton('fs_scale_csv','Download scale results'),
    shiny::h4('Index definitions and limitations'),shiny::verbatimTextOutput('fs_notes'))))
}
ptq_structure_server <- function(input,output,session) {
 for(nm in c('traits','abundance','groups'))local({
  key <- nm
  output[[paste0('fs_',key,'_example')]] <- shiny::downloadHandler(
   filename=function()paste0('functional_structure_',key,'.csv'),
   content=function(file)file.copy(file.path('data',paste0('functional_structure_',key,'.csv')),file,overwrite=TRUE))
 })
 read_csv <- function(upload) {
  shiny::req(upload)
  x <- tryCatch(utils::read.csv(upload$datapath,check.names=FALSE,
   stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
   error=function(e)shiny::validate(shiny::need(FALSE,conditionMessage(e))))
  names(x) <- trimws(names(x))
  shiny::validate(shiny::need(!anyDuplicated(names(x)),'Duplicate CSV columns.'))
  x
 }
 dat <- shiny::reactive({
  t <- read_csv(input$fs_traits_file);a <- read_csv(input$fs_abundance_file)
  shiny::validate(shiny::need('species_id'%in%names(t),'Trait CSV needs species_id.'),
   shiny::need('site_id'%in%names(a),'Abundance CSV needs site_id.'))
  t$species_id <- trimws(as.character(t$species_id))
  a$site_id <- trimws(as.character(a$site_id))
  shiny::validate(shiny::need(!anyNA(t$species_id)&&all(nzchar(t$species_id))&&
   !anyDuplicated(t$species_id),'Invalid species IDs.'),
   shiny::need(!anyNA(a$site_id)&&all(nzchar(a$site_id))&&
   !anyDuplicated(a$site_id),'Invalid site IDs.'))
  species <- setdiff(names(a),'site_id')
  shiny::validate(shiny::need(length(species)>=3,'At least three abundance species columns required.'),
   shiny::need(nrow(a)>=2,'At least two sites required.'))
  mat <- sapply(a[species],function(x)suppressWarnings(as.numeric(as.character(x))))
  mat <- matrix(mat,nrow=nrow(a),dimnames=list(a$site_id,species))
  shiny::validate(shiny::need(all(is.finite(mat))&&all(mat>=0),
   'Abundance values must be nonnegative finite numbers without missing values.'))
  active <- colSums(mat)>0
  mat <- mat[,active,drop=FALSE]
  shiny::validate(shiny::need(ncol(mat)>=3,'At least three species must occur.'),
   shiny::need(all(colnames(mat)%in%t$species_id),
   'Some occurring species have no trait record.'))
  t <- t[match(colnames(mat),t$species_id),,drop=FALSE]
  numeric_traits <- setdiff(names(t),'species_id')
  numeric_traits <- numeric_traits[vapply(t[numeric_traits],function(x){
   v <- suppressWarnings(as.numeric(as.character(x)))
   all(is.na(x)|is.finite(v))
  },logical(1))]
  groups <- NULL
  if(!is.null(input$fs_groups_file)) {
   g <- read_csv(input$fs_groups_file)
   shiny::validate(shiny::need(all(c('site_id','group')%in%names(g)),
    'Site group CSV needs site_id and group.'))
   g$site_id <- trimws(as.character(g$site_id))
   shiny::validate(shiny::need(!anyDuplicated(g$site_id),'Duplicate group site IDs.'),
    shiny::need(all(rownames(mat)%in%g$site_id),'Group CSV missing sites.'))
   groups <- trimws(as.character(g$group[match(rownames(mat),g$site_id)]))
   shiny::validate(shiny::need(all(!is.na(groups)&nzchar(groups)),
    'Missing group labels.'))
  }
  list(t=t,a=mat,traits=numeric_traits,groups=groups)
 })
 output$fs_trait_choices <- shiny::renderUI({
  x <- dat()
  shiny::selectInput('fs_traits','Select numerical traits',choices=x$traits,
   selected=x$traits,multiple=TRUE)
 })
 output$fs_checks <- shiny::renderText({
  if(is.null(input$fs_traits_file)||is.null(input$fs_abundance_file))
   return('Upload paired files or download the examples.')
  x <- dat()
  paste0('Sites: ',nrow(x$a),' | Occurring species: ',ncol(x$a),
   '\nNumeric traits available: ',length(x$traits),
   '\nSite grouping: ',if(is.null(x$groups))'none' else paste(length(unique(x$groups)),'groups'),
   '\nZero-total sites: ',sum(rowSums(x$a)==0))
 })
 calc <- shiny::eventReactive(input$fs_run,{
  x <- dat()
  sel <- intersect(input$fs_traits,x$traits)
  shiny::validate(shiny::need(length(sel)>=2,'Select at least two numeric traits.'))
  T <- sapply(x$t[sel],function(z)as.numeric(as.character(z)))
  T <- matrix(T,nrow=ncol(x$a),dimnames=list(colnames(x$a),sel))
  shiny::validate(shiny::need(all(is.finite(T)),'Missing selected trait values.'),
   shiny::need(all(apply(T,2,stats::sd)>0),'Remove constant traits.'))
  Z <- scale(T)
  D <- as.matrix(stats::dist(Z))
  nsp <- nrow(Z);ns <- nrow(x$a)
  # Species originality: nearest other species distance; specialization:
  # distance to global unweighted trait centroid.
  originality <- vapply(seq_len(nsp),function(i)min(D[i,-i]),numeric(1))
  specialization <- sqrt(rowSums((Z-matrix(colMeans(Z),nsp,ncol(Z),byrow=TRUE))^2))
  sp <- data.frame(species_id=rownames(Z),
   originality_nearest_neighbor=originality,
   specialization_centroid_distance=specialization)
  # Relative abundance based FDis, RaoQ and functional redundancy.
  measure <- function(v) {
   idx <- which(v>0)
   if(!length(idx))return(c(richness=0,FDis=NA_real_,RaoQ=NA_real_,
    redundancy=NA_real_,vulnerability=NA_real_,over_redundancy=NA_real_,
    hull_area_2D=NA_real_))
   p <- v[idx]/sum(v[idx])
   coord <- Z[idx,,drop=FALSE]
   centroid <- colSums(sweep(coord,1,p,'*'))
   fdis <- sum(p*sqrt(rowSums(sweep(coord,2,centroid,'-')^2)))
   rao <- as.numeric(t(p)%*%D[idx,idx,drop=FALSE]%*%p)
   # Redundancy proxy: 1 - RaoQ/max possible species distance,
   # 0..1 for fixed species pool, NOT an independent established index.
   mx <- max(D)
   redundancy <- if(mx>0)1-rao/mx else NA_real_
   # Vulnerability: share of occupied species with no close neighbor,
   # using first quartile of nonzero global interspecific distances.
   threshold <- as.numeric(stats::quantile(D[upper.tri(D)],.25))
   neighbor <- vapply(idx,function(i){
    other <- setdiff(idx,i)
    if(!length(other))0L else sum(D[i,other]<=threshold)
   },integer(1))
   vulnerability <- mean(neighbor==0)
   over <- mean(neighbor>=2)
   hull <- NA_real_
   if(length(idx)>=3) {
    pc <- stats::prcomp(Z,center=FALSE,scale.=FALSE)$x[idx,1:2,drop=FALSE]
    if(nrow(unique(pc))>=3) {
     h <- grDevices::chull(pc[,1],pc[,2]);xx <- pc[h,1];yy <- pc[h,2]
     hull <- abs(sum(xx*c(yy[-1],yy[1])-yy*c(xx[-1],xx[1])))/2
    }
   }
   c(richness=length(idx),FDis=fdis,RaoQ=rao,redundancy=redundancy,
    vulnerability=vulnerability,over_redundancy=over,hull_area_2D=hull)
  }
  obs <- t(apply(x$a,1,measure))
  reps <- as.integer(input$fs_reps);seed <- as.integer(input$fs_seed)
  shiny::validate(shiny::need(is.finite(reps)&&reps>=19&&reps<=999,
   'Permutations must be 19 to 999.'),
   shiny::need(is.finite(seed)&&seed>0,'Seed must be positive.'))
  # Null: shuffle trait identities across species, preserve abundances and occupancy.
  set.seed(seed)
  null_fdis <- matrix(NA_real_,ns,reps)
  null_rao <- matrix(NA_real_,ns,reps)
  for(b in seq_len(reps)) {
   perm <- sample.int(nsp)
   zp <- Z[perm,,drop=FALSE]
   dp <- D[perm,perm,drop=FALSE]
   for(i in seq_len(ns)) {
    v <- x$a[i,];idx <- which(v>0)
    if(!length(idx))next
    p <- v[idx]/sum(v[idx]);zz <- zp[idx,,drop=FALSE]
    c0 <- colSums(sweep(zz,1,p,'*'))
    null_fdis[i,b] <- sum(p*sqrt(rowSums(sweep(zz,2,c0,'-')^2)))
    null_rao[i,b] <- as.numeric(t(p)%*%dp[idx,idx,drop=FALSE]%*%p)
   }
  }
  ses <- function(observed,null) {
   mu <- rowMeans(null,na.rm=TRUE)
   sd <- apply(null,1,stats::sd,na.rm=TRUE)
   z <- (observed-mu)/sd
   z[!is.finite(z)] <- NA_real_
   z
  }
  site <- data.frame(site_id=rownames(x$a),obs,check.names=FALSE)
  site$SES_FDis <- ses(site$FDis,null_fdis)
  site$SES_RaoQ <- ses(site$RaoQ,null_rao)
  site$p_FDis_two_sided <- vapply(seq_len(ns),function(i){
   if(!is.finite(site$FDis[i]))return(NA_real_)
   (1+sum(abs(null_fdis[i,]-mean(null_fdis[i,],na.rm=TRUE))>=
    abs(site$FDis[i]-mean(null_fdis[i,],na.rm=TRUE)),na.rm=TRUE))/(reps+1)
  },numeric(1))
  site$p_RaoQ_two_sided <- vapply(seq_len(ns),function(i){
   if(!is.finite(site$RaoQ[i]))return(NA_real_)
   (1+sum(abs(null_rao[i,]-mean(null_rao[i,],na.rm=TRUE))>=
    abs(site$RaoQ[i]-mean(null_rao[i,],na.rm=TRUE)),na.rm=TRUE))/(reps+1)
  },numeric(1))
  if(!is.null(x$groups))site$group <- x$groups
  # Site accumulation curve: mean pooled richness and FDis over sampled
  # site subsets, up to 50 permutations for responsiveness.
  accum <- list();set.seed(seed+1L)
  for(k in seq_len(ns)) {
   draws <- replicate(min(50L,reps),{
    ids <- sample.int(ns,k)
    m <- measure(colSums(x$a[ids,,drop=FALSE]))
    c(richness=m['richness'],FDis=m['FDis'])
   })
   accum[[k]] <- data.frame(type='accumulation',unit=paste0(k,' sites'),
    n_sites=k,richness=mean(draws[1,],na.rm=TRUE),
    FDis=mean(draws[2,],na.rm=TRUE))
  }
  if(!is.null(x$groups))for(g in unique(x$groups)) {
   ids <- which(x$groups==g)
   m <- measure(colSums(x$a[ids,,drop=FALSE]))
   accum[[length(accum)+1L]] <- data.frame(type='group',unit=g,
    n_sites=length(ids),richness=unname(m['richness']),FDis=unname(m['FDis']))
  }
  list(site=site,species=sp,scale=do.call(rbind,accum),reps=reps)
 })
 output$fs_site <- shiny::renderTable(calc()$site,striped=TRUE,digits=4)
 output$fs_species <- shiny::renderTable(calc()$species,striped=TRUE,digits=4)
 output$fs_scale <- shiny::renderTable(calc()$scale,striped=TRUE,digits=4)
 output$fs_notes <- shiny::renderText({
  paste0('Originality: nearest-neighbor distance in standardized numeric trait space.\n',
   'Specialization: distance from global unweighted species trait centroid.\n',
   'Redundancy: exploratory 1 - RaoQ / maximum global pairwise trait distance, NOT an established functional redundancy estimator.\n',
   'Vulnerability: fraction of present species without neighbors within the global first-quartile distance threshold.\n',
   'Over-redundancy: fraction with at least two such neighbors. These are threshold-based proxies.\n',
   '2D hull area: convex hull on the first two PCA axes; NOT a probabilistic functional hypervolume.\n',
   'SES: (observed - null mean)/null SD, with trait identities shuffled among species; abundance and occupancy fixed.\n',
   'Permutation p-values are exploratory and unadjusted. SES is undefined if null SD=0.\n',
   'Spatial-scale analysis: pooled site-accumulation and optional named group summaries, NOT geospatial multiscale modelling.\n',
   'Trait standardization uses all occurring species; results depend on selected traits and regional species pool.\n',
   'No multidimensional kernel hypervolume, TPD, or formal trait-based redundancy benchmark is implemented.')
 })
 for(name in c('site','species','scale'))local({
  nm <- name
  output[[paste0('fs_',nm,'_csv')]] <- shiny::downloadHandler(
   filename=function()paste0('PlantTraitQuant_functional_',nm,'.csv'),
   content=function(file)utils::write.csv(calc()[[nm]],file,row.names=FALSE,na=''))
 })
}
