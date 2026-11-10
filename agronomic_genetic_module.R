# PlantTraitQuant v1.2.0 | Agronomic trial and genetic variability module
# Native R statistical calculations. Genetic parameters require balanced replication.
ptq_agronomic_ui <- function() {
 shiny::tabPanel('Agronomic & Genetic Variability',
  shiny::h3('Agronomic trait analysis, ANOVA and genetic variability'),
  shiny::p('Upload a long-format trial CSV. Each row is a plot observation. Genotype, replication, environment and numeric traits must be identified.'),
  shiny::fluidRow(
   shiny::column(4,shiny::div(class='well',
    shiny::fileInput('ag_file','Trial observations CSV',accept='.csv'),
    shiny::downloadButton('ag_example','Download example RCBD trial'),
    shiny::downloadButton('ag_multi_example','Download example multi-environment trial'),
    shiny::uiOutput('ag_fields'),
    shiny::selectInput('ag_design','Experimental design',choices=c('Randomized complete block design (RCBD)'='rcbd','Completely randomized design (CRD)'='crd')),
    shiny::selectInput('ag_basis','Variance-component analysis',choices=c('Single environment (balanced RCBD/CRD)'='single','Multi-environment balanced trial'='multi')),
    shiny::numericInput('ag_selection','Selection intensity (%)',5,min=0.1,max=50,step=.1),
    shiny::actionButton('ag_run','Analyse agronomic trial',class='btn-primary'))),
   shiny::column(8,
    shiny::h4('Input and design checks'),shiny::verbatimTextOutput('ag_check'),
    shiny::h4('Descriptive statistics by trait'),shiny::tableOutput('ag_desc'),
    shiny::downloadButton('ag_desc_csv','Download descriptive statistics'),
    shiny::h4('Analysis of variance'),shiny::tableOutput('ag_anova'),
    shiny::downloadButton('ag_anova_csv','Download ANOVA'),
    shiny::h4('Genetic variability and heritability'),shiny::tableOutput('ag_genetics'),
    shiny::downloadButton('ag_genetics_csv','Download genetic parameters'),
    shiny::h4('Genotype mean performance'),shiny::tableOutput('ag_means'),
    shiny::downloadButton('ag_means_csv','Download genotype means'),
    shiny::h4('Trait correlation matrix'),shiny::tableOutput('ag_corr'),
    shiny::downloadButton('ag_corr_csv','Download correlations'),
    shiny::h4('Genotype means plot'),shiny::plotOutput('ag_plot',height='350px'),
    shiny::h4('Interpretation and safeguards'),shiny::verbatimTextOutput('ag_notes'))))
}
ptq_agronomic_server <- function(input,output,session) {
 output$ag_example <- shiny::downloadHandler(
  filename=function()'agronomic_rcbd_example.csv',
  content=function(file)file.copy(file.path('data','agronomic_rcbd_example.csv'),file,overwrite=TRUE))
 output$ag_multi_example <- shiny::downloadHandler(
  filename=function()'agronomic_multi_environment_example.csv',
  content=function(file)file.copy(file.path('data','agronomic_multi_environment_example.csv'),file,overwrite=TRUE))
 dat <- shiny::reactive({
  shiny::req(input$ag_file)
  x <- tryCatch(utils::read.csv(input$ag_file$datapath,check.names=FALSE,
    stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
    error=function(e)shiny::validate(shiny::need(FALSE,conditionMessage(e))))
  names(x) <- trimws(names(x))
  shiny::validate(shiny::need(!anyDuplicated(names(x)),'Column names must be unique.'),
   shiny::need(nrow(x)>=6,'At least six observations required.'))
  x
 })
 output$ag_fields <- shiny::renderUI({
  x <- dat();n <- names(x)
  numeric_candidates <- n[vapply(x,function(z){
   zz <- suppressWarnings(as.numeric(as.character(z)))
   sum(is.finite(zz))>=3
  },logical(1))]
  shiny::tagList(
   shiny::selectInput('ag_genotype','Genotype column',choices=n,selected=if('genotype'%in%n)'genotype' else n[1]),
   shiny::selectInput('ag_rep','Replication / block column',choices=n,selected=if('replication'%in%n)'replication' else n[1]),
   shiny::selectInput('ag_env','Environment (location-year) column',choices=c('None (one environment)'='__single__',stats::setNames(n,n)),selected=if('environment'%in%n)'environment' else '__single__'),
   shiny::selectInput('ag_traits','Numeric agronomic traits',choices=numeric_candidates,selected=intersect(c('plant_height_cm','yield_t_ha','disease_severity_pct'),numeric_candidates),multiple=TRUE)
  )
 })
 prepare <- shiny::reactive({
  x <- dat()
  shiny::req(input$ag_genotype,input$ag_rep,input$ag_traits)
  shiny::validate(shiny::need(length(input$ag_traits)>0,'Select one or more numerical traits.'),
   shiny::need(input$ag_genotype!=input$ag_rep,'Genotype and replication must be different columns.'),
   shiny::need(!any(c(input$ag_genotype,input$ag_rep)%in%input$ag_traits),'Traits must not be ID columns.'))
  x$.gen <- trimws(as.character(x[[input$ag_genotype]]))
  x$.rep <- trimws(as.character(x[[input$ag_rep]]))
  x$.env <- if(identical(input$ag_env,'__single__'))rep('Trial_1',nrow(x)) else trimws(as.character(x[[input$ag_env]]))
  shiny::validate(shiny::need(all(!is.na(x$.gen)&nzchar(x$.gen)),'Genotype IDs cannot be missing.'),
   shiny::need(all(!is.na(x$.rep)&nzchar(x$.rep)),'Replication IDs cannot be missing.'),
   shiny::need(all(!is.na(x$.env)&nzchar(x$.env)),'Environment IDs cannot be missing.'))
  for(tr in input$ag_traits) {
   x[[tr]] <- suppressWarnings(as.numeric(as.character(x[[tr]])))
   shiny::validate(shiny::need(all(is.finite(x[[tr]])),paste('Trait',tr,'contains missing or nonnumeric observations.')))
  }
  key <- paste(x$.env,x$.rep,x$.gen,sep='::')
  shiny::validate(shiny::need(!anyDuplicated(key),'Each environment × replication × genotype combination must occur once. Aggregate subsamples first.'))
  x
 })
 design_check <- shiny::reactive({
  x <- prepare();tab <- table(x$.env,x$.gen)
  rep_tab <- table(x$.env,x$.rep,x$.gen)
  nenv <- length(unique(x$.env));ng <- length(unique(x$.gen))
  # Balanced: every genotype observed in every block within every environment.
  reps_by_env <- tapply(x$.rep,x$.env,function(z)length(unique(z)))
  balanced <- all(rep_tab==1L) && length(unique(as.integer(reps_by_env)))==1L
  if(input$ag_design=='crd') {
   # CRD requires equal numbers of independent replicate units per genotype/env,
   # but does not require each genotype to appear in every named block.
   balanced <- all(tab==tab[1])
  }
  list(nenv=nenv,ng=ng,balanced=balanced,r=as.integer(reps_by_env[1]),
   nobs=nrow(x),envs=unique(x$.env))
 })
 output$ag_check <- shiny::renderText({
  if(is.null(input$ag_file))return('Upload a CSV or download one of the two examples.')
  z <- design_check()
  paste0('Observations: ',z$nobs,' | Genotypes: ',z$ng,' | Environments: ',z$nenv,
   '\nDesign: ',toupper(input$ag_design),' | Replications per first environment: ',z$r,
   '\nBalanced according to selected design: ',z$balanced,
   '\nVariance-component mode: ',input$ag_basis,
   '\nWarning: ANOVA and heritability are design-dependent. Review the assumptions below.')
 })
 calculate <- shiny::eventReactive(input$ag_run,{
  x <- prepare();z <- design_check()
  shiny::validate(shiny::need(z$ng>=2,'At least two genotypes required.'),
   shiny::need(z$r>=2,'At least two replicates required.'),
   shiny::need(z$balanced,'Balanced complete design required for this release. Please verify replication and missing plots.'),
   shiny::need((input$ag_basis=='single'&&z$nenv==1)||(input$ag_basis=='multi'&&z$nenv>=2),
    'Choose the single- or multi-environment variance model to match the uploaded data.'))
  x$.gen <- factor(x$.gen);x$.rep <- factor(x$.rep);x$.env <- factor(x$.env)
  # Enforce same replication count per environment for multi-environment EMS.
  shiny::validate(shiny::need(z$nenv==1||length(unique(as.integer(tapply(as.character(x$.rep),x$.env,function(v)length(unique(v))))))==1,
   'Multi-environment variance components require equal replication counts.'))
  desc <- list();anovas <- list();genetics <- list();means <- list()
  alpha <- as.numeric(input$ag_selection)/100
  # Selection intensity i = phi(z_alpha) / alpha for truncation selection.
  i_sel <- stats::dnorm(stats::qnorm(1-alpha))/alpha
  for(tr in input$ag_traits){
   y <- x[[tr]];avg <- mean(y);sdv <- stats::sd(y)
   desc[[tr]] <- data.frame(trait=tr,n=length(y),mean=avg,sd=sdv,
    variance=stats::var(y),minimum=min(y),maximum=max(y),range=diff(range(y)),
    CV_percent=if(abs(avg)>1e-12)100*sdv/abs(avg) else NA_real_)
   dd <- data.frame(y=y,genotype=x$.gen,replication=x$.rep,environment=x$.env)
   if(z$nenv==1) {
    if(input$ag_design=='rcbd')fit <- stats::lm(y~replication+genotype,data=dd)
    else fit <- stats::lm(y~genotype,data=dd)
   } else {
    if(input$ag_design=='rcbd')fit <- stats::lm(y~environment+environment:replication+genotype+environment:genotype,data=dd)
    else fit <- stats::lm(y~environment+genotype+environment:genotype,data=dd)
   }
   at <- stats::anova(fit)
   at$term <- rownames(at)
   at$trait <- tr
   anovas[[tr]] <- data.frame(trait=tr,term=at$term,df=at$Df,
    sum_sq=at[['Sum Sq']],mean_sq=at[['Mean Sq']],F_value=at[['F value']],p_value=at[['Pr(>F)']])
   getms <- function(pattern){
    j <- which(rownames(at)==pattern)
    if(length(j)==1)as.numeric(at[['Mean Sq']][j]) else NA_real_
   }
   mserr <- getms('Residuals');msg <- getms('genotype')
   if(z$nenv==1){
    vg_raw <- (msg-mserr)/z$r
    vge_raw <- NA_real_
    vp <- vg_raw+mserr
    h2_denom <- vg_raw+mserr/z$r
   } else {
    msge <- getms('environment:genotype')
    vg_raw <- (msg-msge)/(z$r*z$nenv)
    vge_raw <- (msge-mserr)/z$r
    vp <- vg_raw+vge_raw+mserr
    h2_denom <- vg_raw+vge_raw/z$nenv+mserr/(z$r*z$nenv)
   }
   # Negative ANOVA method-of-moments estimates are truncated for derived measures,
   # but retained in raw columns to make boundary handling transparent.
   vg <- max(0,vg_raw);vge <- if(is.na(vge_raw))NA_real_ else max(0,vge_raw)
   if(z$nenv==1)vp <- vg+mserr
   else vp <- vg+vge+mserr
   h2den <- if(z$nenv==1)vg+mserr/z$r else vg+vge/z$nenv+mserr/(z$r*z$nenv)
   h2 <- if(h2den>0)vg/h2den else NA_real_
   gcv <- if(abs(avg)>1e-12)100*sqrt(vg)/abs(avg) else NA_real_
   pcv <- if(abs(avg)>1e-12)100*sqrt(vp)/abs(avg) else NA_real_
   ga <- if(is.finite(h2))i_sel*sqrt(h2den)*h2 else NA_real_
   genetics[[tr]] <- data.frame(trait=tr,genotypic_variance_raw=vg_raw,
    genotypic_variance_truncated=vg,genotype_environment_variance_raw=vge_raw,
    residual_variance=mserr,phenotypic_plot_variance=vp,
    GCV_percent=gcv,PCV_percent=pcv,
    H2_entry_mean=h2,H2_percent=100*h2,
    genetic_advance=ga,GA_percent_mean=if(abs(avg)>1e-12)100*ga/abs(avg) else NA_real_,
    selection_intensity_percent=100*alpha)
   tmp <- stats::aggregate(dd$y,by=list(genotype=dd$genotype),FUN=mean)
   names(tmp)[2] <- 'genotype_mean';tmp$trait <- tr
   means[[tr]] <- tmp[,c('trait','genotype','genotype_mean')]
  }
  # Correlations computed on genotype means, not individual plot observations.
  mm <- stats::aggregate(x[input$ag_traits],by=list(genotype=x$.gen),FUN=mean)
  cm <- if(length(input$ag_traits)>1)stats::cor(mm[input$ag_traits]) else
   matrix(1,1,1,dimnames=list(input$ag_traits,input$ag_traits))
  corr <- data.frame(trait=rownames(cm),cm,check.names=FALSE,row.names=NULL)
  list(desc=do.call(rbind,desc),anova=do.call(rbind,anovas),
   genetics=do.call(rbind,genetics),means=do.call(rbind,means),corr=corr,
   mode=input$ag_basis,design=input$ag_design,nenv=z$nenv,r=z$r)
 })
 for(nm in c('desc','anova','genetics','means','corr'))local({
  key <- nm
  output[[paste0('ag_',key)]] <- shiny::renderTable(calculate()[[key]],striped=TRUE,digits=4)
  output[[paste0('ag_',key,'_csv')]] <- shiny::downloadHandler(
   filename=function()paste0('PlantTraitQuant_agronomic_',key,'.csv'),
   content=function(file)utils::write.csv(calculate()[[key]],file,row.names=FALSE,na=''))
 })
 output$ag_plot <- shiny::renderPlot({
  a <- calculate()$means
  a <- a[a$trait==a$trait[1],,drop=FALSE]
  o <- order(a$genotype_mean,decreasing=TRUE)
  graphics::barplot(a$genotype_mean[o],names.arg=as.character(a$genotype[o]),
   las=2,cex.names=.7,main=paste('Genotype means:',a$trait[1]),
   ylab=a$trait[1],col='#2c7a57',border=NA)
 })
 output$ag_notes <- shiny::renderText({
  z <- calculate()
  paste0('Model: ',if(z$nenv==1)'single environment' else 'balanced multi-environment',
   '; design: ',toupper(z$design),'; replications: ',z$r,'.\n',
   'Genotypic variance from ANOVA expected mean squares; negative raw estimates are truncated to zero only for derived indices.\n',
   'Single environment: sigma_g^2=(MS_genotype-MS_error)/r.\n',
   'Multi-environment: sigma_g^2=(MS_genotype-MS_GxE)/(r*e); sigma_GxE^2=(MS_GxE-MS_error)/r.\n',
   'GCV=100*sqrt(sigma_g^2)/abs(grand mean); PCV uses plot-level sigma_g^2 + sigma_GxE^2 + sigma_error^2.\n',
   'Broad-sense H2 is on an entry-mean basis: sigma_g^2/[sigma_g^2+sigma_GxE^2/e+sigma_error^2/(r*e)].\n',
   'Genetic advance = selection intensity * sqrt(entry-mean phenotypic variance) * H2.\n',
   'Variance components are method-of-moments approximations for balanced experiments, not REML estimates.\n',
   'Genotype effects are interpreted as random for variance components, although ANOVA tests use fixed-effect lm; this distinction requires care.\n',
   'For multi-environment RCBD, replication is nested within environment; for CRD, replication is an independent experimental unit.\n',
   'No estimates of narrow-sense heritability, BLUPs, genomic heritability, or adjusted genotype means are claimed.\n',
   'Correlation matrix is calculated from unadjusted genotype means.\n',
   'Check independence, randomization, homoscedasticity, residual normality and environmental comparability before inference.')
 })
}
