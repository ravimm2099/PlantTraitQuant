# Full ade4 RLQ and fourth-corner interface. Optional dependency: ade4.
ptq_ade4_ui <- function() {
 shiny::tabPanel('RLQ & Fourth-Corner',
  shiny::h3('RLQ ordination and fourth-corner tests'),
  shiny::p('Reference-package workflow: correspondence analysis of L, weighted PCA of R and Q, ade4::rlq, and ade4::fourthcorner.'),
  shiny::fluidRow(
   shiny::column(4,shiny::div(class='well',
    shiny::fileInput('a4_l','L: Site-by-species abundance',accept='.csv'),
    shiny::fileInput('a4_r','R: Site environmental variables',accept='.csv'),
    shiny::fileInput('a4_q','Q: Species functional traits',accept='.csv'),
    shiny::numericInput('a4_perm','Permutations for fourth-corner tests',199,min=19,max=9999),
    shiny::numericInput('a4_seed','Random seed',123,min=1,max=1000000),
    shiny::actionButton('a4_run','Run ade4 analysis',class='btn-primary'),
    shiny::downloadButton('a4_example_l','Example L'),
    shiny::downloadButton('a4_example_r','Example R'),
    shiny::downloadButton('a4_example_q','Example Q'),
    shiny::p('Requires optional package ade4. Install using install.packages("ade4").'))),
   shiny::column(8,
    shiny::h4('Validation and package status'),shiny::verbatimTextOutput('a4_checks'),
    shiny::h4('RLQ eigenvalues'),shiny::tableOutput('a4_eigen'),
    shiny::h4('RLQ environmental loadings'),shiny::tableOutput('a4_rload'),
    shiny::h4('RLQ trait loadings'),shiny::tableOutput('a4_qload'),
    shiny::plotOutput('a4_plot',height='420px'),
    shiny::downloadButton('a4_eigen_csv','Download eigenvalues'),
    shiny::downloadButton('a4_rload_csv','Download R loadings'),
    shiny::downloadButton('a4_qload_csv','Download Q loadings'),
    shiny::h4('Fourth-corner reference-package output'),
    shiny::verbatimTextOutput('a4_fc'),
    shiny::downloadButton('a4_fc_txt','Download fourth-corner output'),
    shiny::h4('Interpretation'),shiny::verbatimTextOutput('a4_notes'))))
}
ptq_ade4_server <- function(input,output,session) {
 for(nm in c('l','r','q'))local({
  name <- nm
  output[[paste0('a4_example_',name)]] <- shiny::downloadHandler(
   filename=function()paste0('trait_environment_',name,'.csv'),
   content=function(file)file.copy(file.path('data',paste0('trait_environment_',name,'.csv')),file,overwrite=TRUE))
 })
 read_one <- function(upload,id) {
  shiny::req(upload)
  d <- tryCatch(utils::read.csv(upload$datapath,check.names=FALSE,
   stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
   error=function(e)shiny::validate(shiny::need(FALSE,conditionMessage(e))))
  names(d) <- trimws(names(d))
  shiny::validate(shiny::need(!anyDuplicated(names(d)),'Duplicate CSV column names.'),
   shiny::need(id%in%names(d),paste('Required column missing:',id)))
  d[[id]] <- trimws(as.character(d[[id]]))
  shiny::validate(shiny::need(!anyNA(d[[id]]) && all(nzchar(d[[id]])) &&
   !anyDuplicated(d[[id]]),paste('Invalid or duplicated',id)))
  d
 }
 tables <- shiny::reactive({
  L <- read_one(input$a4_l,'site_id')
  R <- read_one(input$a4_r,'site_id')
  Q <- read_one(input$a4_q,'species_id')
  species <- setdiff(names(L),'site_id')
  shiny::validate(shiny::need(length(species)>=2,'L requires at least two species.'),
   shiny::need(all(L$site_id%in%R$site_id),'R is missing site IDs present in L.'),
   shiny::need(all(species%in%Q$species_id),'Q is missing species IDs present in L.'))
  R <- R[match(L$site_id,R$site_id),,drop=FALSE]
  Q <- Q[match(species,Q$species_id),,drop=FALSE]
  numeric_df <- function(d,id) {
   cols <- setdiff(names(d),id)
   shiny::validate(shiny::need(length(cols)>=2,paste(id,'table requires at least two numeric variables.')))
   v <- lapply(d[cols],function(x)suppressWarnings(as.numeric(as.character(x))))
   z <- as.data.frame(v,check.names=FALSE)
   names(z) <- cols
   shiny::validate(shiny::need(all(vapply(z,function(x)all(is.finite(x)),logical(1))),
    paste(id,'table has missing/non-numeric data. Complete numeric data are required.')),
    shiny::need(all(vapply(z,function(x)stats::sd(x)>0,logical(1))),
     paste(id,'table contains constant variables.')))
   z
  }
  a <- as.matrix(numeric_df(L,'site_id'))
  shiny::validate(shiny::need(all(a>=0),'L abundances must be nonnegative.'),
   shiny::need(nrow(a)>=4,'At least four sites are required.'),
   shiny::need(all(rowSums(a)>0) && all(colSums(a)>0),
    'L must have positive abundance in every site and for every species.'))
  rr <- numeric_df(R,'site_id');qq <- numeric_df(Q,'species_id')
  rownames(rr) <- L$site_id;rownames(qq) <- species
  rownames(a) <- L$site_id;colnames(a) <- species
  list(L=as.data.frame(a,check.names=FALSE),R=rr,Q=qq)
 })
 output$a4_checks <- shiny::renderText({
  if(!requireNamespace('ade4',quietly=TRUE))
   return('ade4 is not installed. Run install.packages("ade4") in RStudio.')
  if(is.null(input$a4_l)||is.null(input$a4_r)||is.null(input$a4_q))
   return(paste('ade4 version:',as.character(utils::packageVersion('ade4')),
    '\nUpload three CSVs or download the examples.'))
  x <- tables()
  paste0('ade4 version: ',as.character(utils::packageVersion('ade4')),
   '\nSites: ',nrow(x$L),' | Species: ',ncol(x$L),
   '\nEnvironment variables: ',ncol(x$R),' | Traits: ',ncol(x$Q),
   '\nRows and columns matched by IDs, not original file order.')
 })
 analysis <- shiny::eventReactive(input$a4_run,{
  shiny::validate(shiny::need(requireNamespace('ade4',quietly=TRUE),
   'Please install ade4: install.packages("ade4")'))
  x <- tables()
  permutations <- as.integer(input$a4_perm)
  seed <- as.integer(input$a4_seed)
  shiny::validate(shiny::need(!is.na(permutations) && permutations>=19 && permutations<=9999,
   'Permutations must be 19 to 9999.'),
   shiny::need(!is.na(seed) && seed>0,'Seed must be positive.'))
  set.seed(seed)
  tryCatch({
   # Standard ade4 three-table approach, with CA row/column weights.
   coa <- ade4::dudi.coa(x$L,scannf=FALSE,nf=2)
   dR <- ade4::dudi.pca(x$R,row.w=coa$lw,center=TRUE,scale=TRUE,scannf=FALSE,nf=2)
   dQ <- ade4::dudi.pca(x$Q,row.w=coa$cw,center=TRUE,scale=TRUE,scannf=FALSE,nf=2)
   rlq <- ade4::rlq(dR,coa,dQ,scannf=FALSE,nf=2)
   # Model 2 permutes sites; model 4 permutes species. They are separate null hypotheses.
   fc2 <- ade4::fourthcorner(x$R,x$L,x$Q,modeltype=2,nrepet=permutations)
   fc4 <- ade4::fourthcorner(x$R,x$L,x$Q,modeltype=4,nrepet=permutations)
   list(rlq=rlq,fc2=fc2,fc4=fc4,coa=coa,
    version=as.character(utils::packageVersion('ade4')))
  },error=function(e)shiny::validate(shiny::need(FALSE,
   paste('ade4 analysis failed:',conditionMessage(e)))))
 })
 eigen <- shiny::reactive({
  z <- analysis()$rlq$eig
  data.frame(axis=seq_along(z),eigenvalue=as.numeric(z),
   explained_percent=100*as.numeric(z)/sum(z))
 })
 rload <- shiny::reactive({
  x <- analysis()$rlq$l1
  if(is.null(x))return(data.frame(message='RLQ R loadings not available'))
  data.frame(variable=rownames(x),x,check.names=FALSE,row.names=NULL)
 })
 qload <- shiny::reactive({
  x <- analysis()$rlq$c1
  if(is.null(x))return(data.frame(message='RLQ Q loadings not available'))
  data.frame(variable=rownames(x),x,check.names=FALSE,row.names=NULL)
 })
 output$a4_eigen <- shiny::renderTable(eigen(),striped=TRUE,digits=4)
 output$a4_rload <- shiny::renderTable(rload(),striped=TRUE,digits=4)
 output$a4_qload <- shiny::renderTable(qload(),striped=TRUE,digits=4)
 output$a4_plot <- shiny::renderPlot({
  z <- analysis()$rlq
  scores <- z$l1
  if(is.null(scores)||ncol(scores)<2) {
   graphics::plot.new();graphics::text(.5,.5,'Two RLQ axes are not available.')
  } else {
   graphics::plot(scores[,1],scores[,2],pch=19,
    xlab='RLQ axis 1',ylab='RLQ axis 2',main='RLQ environmental variable scores')
   graphics::text(scores[,1],scores[,2],labels=rownames(scores),pos=3,cex=.8)
   graphics::abline(h=0,v=0,col='grey80')
  }
 })
 fc_text <- shiny::reactive({
  z <- analysis()
  paste(c(paste('ade4',z$version,'| Model 2: permute sites'),
   utils::capture.output(print(z$fc2)),
   '\nModel 4: permute species',utils::capture.output(print(z$fc4))),collapse='\n')
 })
 output$a4_fc <- shiny::renderText(fc_text())
 output$a4_notes <- shiny::renderText({
  shiny::req(analysis())
  paste('This tab calls ade4::rlq and ade4::fourthcorner directly.',
   'Model 2 and Model 4 outputs are reported separately, not combined into a model-6 test.',
   'Permutation p-values depend on the exchangeability assumptions and spatial structure.',
   'The RLQ loading tables are native ade4 outputs, not independently recomputed estimates.',
   'This version has not been numerically benchmarked against a separate reference dataset.',sep='\n')
 })
 output$a4_eigen_csv <- shiny::downloadHandler(
  filename=function()'PlantTraitQuant_ade4_RLQ_eigenvalues.csv',
  content=function(file)utils::write.csv(eigen(),file,row.names=FALSE))
 output$a4_rload_csv <- shiny::downloadHandler(
  filename=function()'PlantTraitQuant_ade4_R_loadings.csv',
  content=function(file)utils::write.csv(rload(),file,row.names=FALSE))
 output$a4_qload_csv <- shiny::downloadHandler(
  filename=function()'PlantTraitQuant_ade4_Q_loadings.csv',
  content=function(file)utils::write.csv(qload(),file,row.names=FALSE))
 output$a4_fc_txt <- shiny::downloadHandler(
  filename=function()'PlantTraitQuant_ade4_fourthcorner.txt',
  content=function(file)writeLines(fc_text(),file))
}
