
# Diploid biallelic SNP dosage workflow (0/1/2 alternate-allele copies).
ptq_snp_ui <- function() {
 shiny::tabPanel('SNP diversity & genomic core',
  shiny::h3('Diploid biallelic SNP diversity and core collection'),
  shiny::p('Input: accession_id and SNP marker columns with 0, 1 or 2 alternate-allele copies. Missing calls may be blank or NA.'),
  shiny::fluidRow(
   shiny::column(4,shiny::div(class='well',
    shiny::fileInput('snp_file','SNP dosage matrix CSV',accept='.csv'),
    shiny::numericInput('snp_maf','Minimum minor allele frequency',0,min=0,max=0.5,step=0.01),
    shiny::numericInput('snp_call','Minimum marker call rate',0.8,min=0,max=1,step=0.05),
    shiny::numericInput('snp_core','Core size (%)',20,min=1,max=100),
    shiny::numericInput('snp_mini','Mini-core size (%)',5,min=1,max=100),
    shiny::numericInput('snp_seed','Tie-breaking seed',123,min=1,max=1000000),
    shiny::actionButton('snp_go','Analyse SNP data',class='btn-primary'),
    shiny::downloadButton('snp_example','Download example SNP CSV'))),
   shiny::column(8,
    shiny::h4('Validation and filtering'),shiny::verbatimTextOutput('snp_checks'),
    shiny::h4('Marker diversity'),shiny::tableOutput('snp_stats'),
    shiny::downloadButton('snp_stats_csv','Download marker statistics'),
    shiny::h4('Genomic core and mini-core coverage'),shiny::tableOutput('snp_summary'),
    shiny::h4('Selected accessions'),shiny::tableOutput('snp_selected'),
    shiny::downloadButton('snp_selected_csv','Download accession selection'),
    shiny::downloadButton('snp_summary_csv','Download retention summary'),
    shiny::h4('Scientific scope'),shiny::verbatimTextOutput('snp_notes'))))
}
ptq_snp_server <- function(input,output,session) {
 output$snp_example <- shiny::downloadHandler(
  filename=function()'diploid_snp_example.csv',
  content=function(file)file.copy(file.path('data','diploid_snp_example.csv'),file,overwrite=TRUE))
 dat <- shiny::reactive({
  shiny::req(input$snp_file)
  x <- tryCatch(utils::read.csv(input$snp_file$datapath,check.names=FALSE,
   stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
   error=function(e)shiny::validate(shiny::need(FALSE,conditionMessage(e))))
  names(x) <- trimws(names(x))
  shiny::validate(shiny::need('accession_id'%in%names(x),'Missing accession_id column.'),
   shiny::need(!anyDuplicated(names(x)),'Duplicate marker names.'),
   shiny::need(nrow(x)>=3 && nrow(x)<=500,'Requires 3 to 500 accessions.'))
  id <- trimws(as.character(x$accession_id))
  shiny::validate(shiny::need(!anyNA(id) && all(nzchar(id)) && !anyDuplicated(id),
   'Accession IDs must be unique and nonempty.'))
  names0 <- setdiff(names(x),'accession_id')
  shiny::validate(shiny::need(length(names0)>=2,'Provide at least two SNP markers.'))
  a <- matrix(NA_real_,nrow=nrow(x),ncol=length(names0),dimnames=list(id,names0))
  for(j in seq_along(names0)) {
   v <- trimws(as.character(x[[names0[j]]]))
   v[v %in% c('','NA','NaN')] <- NA_character_
   parsed <- suppressWarnings(as.numeric(v))
   shiny::validate(shiny::need(all(is.na(v) | (!is.na(parsed) & parsed%in%c(0,1,2))),
    paste('Invalid genotype in',names0[j],': only 0, 1, 2, or missing are allowed.')))
   a[,j] <- parsed
  }
  a
 })
 output$snp_checks <- shiny::renderText({
  if(is.null(input$snp_file))return('Upload a SNP dosage matrix or download the example.')
  a <- dat()
  paste0('Accessions: ',nrow(a),' | Input SNP markers: ',ncol(a),
    '\nMissing genotype cells: ',sum(is.na(a)),
    '\nGenotype codes: 0, 1, 2 (alternate allele dosage).',
    '\nFiltering applied after Analyse SNP data.')
 })
 res <- shiny::eventReactive(input$snp_go,{
  a <- dat();n <- nrow(a)
  maf_min <- input$snp_maf;call_min <- input$snp_call
  shiny::validate(shiny::need(is.finite(maf_min) && maf_min>=0 && maf_min<=0.5,
   'MAF threshold must be 0 to 0.5.'),
   shiny::need(is.finite(call_min) && call_min>=0 && call_min<=1,
   'Call rate threshold must be 0 to 1.'),
   shiny::need(is.finite(input$snp_core) && input$snp_core>0 && input$snp_core<=100,
   'Core percentage must be 1 to 100.'),
   shiny::need(is.finite(input$snp_mini) && input$snp_mini>0 &&
    input$snp_mini<=input$snp_core,'Mini-core must be positive and no larger than core.'))
  callrate <- colMeans(!is.na(a))
  alt <- colMeans(a,na.rm=TRUE)/2
  alt[!is.finite(alt)] <- NA_real_
  maf <- pmin(alt,1-alt)
  keep <- is.finite(maf) & maf>0 & maf>=maf_min & callrate>=call_min
  shiny::validate(shiny::need(any(keep),
   'No polymorphic markers pass call-rate and MAF filtering.'))
  b <- a[,keep,drop=FALSE]; p <- alt[keep];m <- ncol(b)
  stats <- data.frame(marker=colnames(b),n_called=colSums(!is.na(b)),
   call_rate=callrate[keep],alt_allele_frequency=p,MAF=maf[keep],
   Ho=colMeans(b==1,na.rm=TRUE),He=2*p*(1-p),
   PIC=2*p*(1-p)*(1-2*p*(1-p)),row.names=NULL)
  # For each polymorphic SNP, retain both observed allelic states.
  # A dosage 1 contains both; 0 contains reference only; 2 alternate only.
  features <- cbind(b<2,b>0)
  features[is.na(features)] <- FALSE
  colnames(features) <- c(paste0(colnames(b),'::REF'),
                          paste0(colnames(b),'::ALT'))
  # Pairwise distance: mean absolute dosage difference / 2 over jointly called SNPs.
  d <- matrix(1,n,n)
  diag(d) <- 0
  for(i in seq_len(n))if(i<n)for(j in (i+1):n) {
   valid <- !is.na(b[i,]) & !is.na(b[j,])
   v <- if(any(valid))mean(abs(b[i,valid]-b[j,valid]))/2 else 1
   d[i,j] <- d[j,i] <- v
  }
  seed <- as.integer(input$snp_seed)
  shiny::validate(shiny::need(is.finite(seed) && seed>=1,'Seed must be positive.'))
  set.seed(seed)
  choose <- function(pool,k) {
   selected <- integer()
   while(length(selected)<k) {
    candidates <- setdiff(pool,selected)
    covered <- if(length(selected))colSums(features[selected,,drop=FALSE])>0 else
     rep(FALSE,ncol(features))
    gain <- vapply(candidates,function(i)sum(features[i,] & !covered),numeric(1))
    candidates <- candidates[gain==max(gain)]
    if(length(candidates)>1 && length(selected)) {
     spread <- vapply(candidates,function(i)min(d[i,selected]),numeric(1))
     candidates <- candidates[spread==max(spread)]
    }
    selected <- c(selected,sample(candidates,1))
   }
   selected
  }
  kc <- min(n,max(1L,ceiling(n*input$snp_core/100)))
  km <- min(kc,max(1L,ceiling(n*input$snp_mini/100)))
  core <- choose(seq_len(n),kc)
  mini <- choose(core,km)
  total <- sum(colSums(features)>0)
  retained <- function(idx)sum(colSums(features[idx,,drop=FALSE])>0)
  summary <- data.frame(collection=c('Core','Mini-core'),
   accessions=c(kc,km),percent_of_accessions=100*c(kc,km)/n,
   observed_allele_states_retained=c(retained(core),retained(mini)),
   observed_allele_states_total=total,
   allele_state_retention_percent=100*c(retained(core),retained(mini))/total)
  selection <- data.frame(accession_id=rownames(a),
   collection=ifelse(seq_len(n)%in%mini,'Mini-core (also core)',
    ifelse(seq_len(n)%in%core,'Core only','Not selected')))
  list(stats=stats,summary=summary,selection=selection,kept=m,
   removed=ncol(a)-m)
 })
 output$snp_stats <- shiny::renderTable(utils::head(res()$stats,100),digits=4,striped=TRUE)
 output$snp_summary <- shiny::renderTable(res()$summary,digits=3,striped=TRUE)
 output$snp_selected <- shiny::renderTable(
  res()$selection[res()$selection$collection!='Not selected',,drop=FALSE],striped=TRUE)
 output$snp_notes <- shiny::renderText({
  z <- res()
  paste0('Retained polymorphic SNPs: ',z$kept,'; excluded markers: ',z$removed,
   '\nHo = fraction of observed heterozygotes. He = 2p(1-p) (uncorrected).',
   '\nPIC = 2p(1-p)[1-2p(1-p)] for biallelic codominant markers.',
   '\nCore selection greedily maximizes observed reference/alternate allele-state coverage.',
   '\nGenotype distance is mean absolute dosage difference / 2 at jointly called SNPs.',
   '\nUncalled genotypes contribute neither allele state; unshared pairs receive distance 1.',
   '\nMini-core is nested within core. No global optimum or independent numerical benchmark is claimed.',
   '\nNo LD pruning, kinship adjustment, population structure, Hardy-Weinberg testing, imputation, or rare-allele weighting.',
   '\nSNP allele-state retention can saturate quickly; it does not establish genomic representativeness.')
 })
 output$snp_stats_csv <- shiny::downloadHandler(
  filename=function()'PlantTraitQuant_SNP_marker_statistics.csv',
  content=function(file)utils::write.csv(res()$stats,file,row.names=FALSE))
 output$snp_selected_csv <- shiny::downloadHandler(
  filename=function()'PlantTraitQuant_SNP_core_selection.csv',
  content=function(file)utils::write.csv(res()$selection,file,row.names=FALSE))
 output$snp_summary_csv <- shiny::downloadHandler(
  filename=function()'PlantTraitQuant_SNP_allele_state_retention.csv',
  content=function(file)utils::write.csv(res()$summary,file,row.names=FALSE))
}
