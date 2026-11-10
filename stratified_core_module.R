# Additional accession selection workflows. Optional published packages are called through namespaces.
# PlantTraitQuant does not rebrand their algorithms as original methods.
ptq_stratified_ui <- function() {
 shiny::tabPanel('Stratified Core Selection',
  shiny::h3('Stratified core and mini-core collection workflows'),
  shiny::p('Upload one accession descriptor CSV. It must contain accession_id, group and one or more descriptors. Use optional published SampleCore or Core Hunter backends when installed.'),
  shiny::fluidRow(
   shiny::column(4,shiny::div(class='well',
    shiny::fileInput('sc_data','Accession descriptor CSV',accept='.csv'),
    shiny::downloadButton('sc_example','Download example descriptors'),
    shiny::uiOutput('sc_columns'),
    shiny::sliderInput('sc_fraction','Core fraction',min=.05,max=.5,value=.2,step=.01),
    shiny::selectInput('sc_allocation','Group allocation',choices=c('Proportional'='prop','Equal across groups'='const','Square root'='sqrt','Logarithmic'='log','Distance-weighted'='distance')),
    shiny::selectInput('sc_method','Selection strategy',choices=c('Stratified random (native)'='random','Stratified farthest-first (native)'='maximin','SampleCore random (published package)'='samplecore','Core Hunter phenotype optimization (published package)'='corehunter')),
    shiny::numericInput('sc_seed','Random seed',123,min=1,max=1000000),
    shiny::actionButton('sc_run','Construct core',class='btn-primary'))),
   shiny::column(8,
    shiny::h4('Validation and package availability'),shiny::verbatimTextOutput('sc_check'),
    shiny::h4('Allocation by group'),shiny::tableOutput('sc_alloc'),
    shiny::downloadButton('sc_alloc_csv','Download allocation'),
    shiny::h4('Selected core accessions'),shiny::tableOutput('sc_selected'),
    shiny::downloadButton('sc_selected_csv','Download selected accessions'),
    shiny::h4('Full vs core evaluation'),shiny::tableOutput('sc_eval'),
    shiny::downloadButton('sc_eval_csv','Download evaluation'),
    shiny::h4('Method notes'),shiny::verbatimTextOutput('sc_notes'))))
}
ptq_stratified_server <- function(input,output,session) {
 output$sc_example <- shiny::downloadHandler(
  filename=function()'stratified_core_example.csv',
  content=function(file)file.copy(file.path('data','stratified_core_example.csv'),file,overwrite=TRUE))
 dat <- shiny::reactive({
  shiny::req(input$sc_data)
  x <- tryCatch(utils::read.csv(input$sc_data$datapath,check.names=FALSE,
   stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
   error=function(e)shiny::validate(shiny::need(FALSE,conditionMessage(e))))
  names(x) <- trimws(names(x))
  shiny::validate(shiny::need(!anyDuplicated(names(x)),'Duplicate column names'),
   shiny::need('accession_id'%in%names(x),'Required column: accession_id'),
   shiny::need('group'%in%names(x),'Required column: group'),
   shiny::need(nrow(x)>=4,'At least four accessions required'))
  x$accession_id <- trimws(as.character(x$accession_id))
  x$group <- trimws(as.character(x$group))
  shiny::validate(shiny::need(!anyNA(x$accession_id)&&all(nzchar(x$accession_id))&&
   !anyDuplicated(x$accession_id),'Unique nonempty accession IDs required'),
   shiny::need(!anyNA(x$group)&&all(nzchar(x$group)),'Group must be nonempty'))
  x
 })
 output$sc_columns <- shiny::renderUI({
  x<-dat();candidates<-setdiff(names(x),c('accession_id','group'))
  shiny::selectInput('sc_traits','Descriptors for selection and evaluation',choices=candidates,
   selected=candidates,multiple=TRUE)
 })
 output$sc_check <- shiny::renderText({
  if(is.null(input$sc_data))return('Upload a descriptor CSV or download the example.')
  x<-dat()
  paste('Accessions:',nrow(x),'| Groups:',length(unique(x$group)),
   '\nSampleCore installed:',requireNamespace('SampleCore',quietly=TRUE),
   '\ncorehunter installed:',requireNamespace('corehunter',quietly=TRUE),
   '\nCore Hunter also requires Java and rJava.')
 })
 # Allocate exactly n entries subject to group capacities using a capped largest-remainder method.
 allocate <- function(groups,n,method,D=NULL) {
  counts<-table(groups);capacity<-as.integer(counts);names(capacity)<-names(counts)
  g<-length(capacity)
  if(method=='const')weight<-rep(1,g)
  else if(method=='sqrt')weight<-sqrt(capacity)
  else if(method=='log')weight<-log1p(capacity)
  else if(method=='distance')weight<-pmax(D[names(capacity)],0)
  else weight<-as.numeric(capacity)
  if(!all(is.finite(weight))||sum(weight)<=0)weight<-as.numeric(capacity)
  allocated<-integer(g);names(allocated)<-names(capacity)
  # Include all groups if size allows; prevents accidental omission of small groups.
  if(n>=g){allocated[]<-1L;nleft<-n-g}else nleft<-n
  while(nleft>0){
   eligible<-which(allocated<capacity)
   if(!length(eligible))break
   score<-weight[eligible]/(allocated[eligible]+1)
   j<-eligible[which.max(score)]
   allocated[j]<-allocated[j]+1L;nleft<-nleft-1L
  }
  allocated
 }
 distmat <- function(x,traits) {
  shiny::validate(shiny::need(length(traits)>0,'Select at least one descriptor'))
  y<-x[,traits,drop=FALSE]
  for(nm in names(y)){
   if(is.character(y[[nm]]))y[[nm]]<-factor(y[[nm]])
   if(is.logical(y[[nm]]))y[[nm]]<-factor(y[[nm]])
  }
  shiny::validate(shiny::need(!anyNA(y),'Selected descriptors cannot contain missing values'),
   shiny::need(requireNamespace('cluster',quietly=TRUE),'Install package cluster for Gower distances'))
  as.matrix(cluster::daisy(y,metric='gower'))
 }
 choose_maximin <- function(idx,k,D) {
  if(k<=0)return(integer(0))
  if(k>=length(idx))return(idx)
  # Deterministic within-group farthest-first seed: max mean distance.
  sub<-D[idx,idx,drop=FALSE]
  chosen<-idx[which.max(rowMeans(sub))]
  while(length(chosen)<k){
   remaining<-setdiff(idx,chosen)
   score<-vapply(remaining,function(i)min(D[i,chosen]),numeric(1))
   chosen<-c(chosen,remaining[which.max(score)])
  }
  chosen
 }
 result <- shiny::eventReactive(input$sc_run,{
  x<-dat();traits<-intersect(input$sc_traits,setdiff(names(x),c('accession_id','group')))
  shiny::validate(shiny::need(length(traits)>0,'Select at least one descriptor'))
  n<-max(1L,min(nrow(x),as.integer(round(nrow(x)*input$sc_fraction))))
  method<-input$sc_method;strategy<-input$sc_allocation
  shiny::validate(shiny::need(is.finite(input$sc_seed),'Set a valid random seed'))
  set.seed(as.integer(input$sc_seed))
  D<-NULL
  if(strategy=='distance'||method=='maximin')D<-distmat(x,traits)
  weights<-NULL
  if(strategy=='distance'){
   groups<-split(seq_len(nrow(x)),x$group)
   weights<-vapply(groups,function(ix){
    if(length(ix)<2)return(0)
    d<-D[ix,ix,drop=FALSE];mean(d[upper.tri(d)])
   },numeric(1))
  }
  alloc<-allocate(x$group,n,strategy,weights)
  allocation<-data.frame(group=names(alloc),total=as.integer(table(x$group)[names(alloc)]),selected=as.integer(alloc))
  if(method=='corehunter'){
   shiny::validate(shiny::need(requireNamespace('corehunter',quietly=TRUE),
    'Install the supplied corehunter package and its Java/rJava dependencies first.'))
   y<-x[,traits,drop=FALSE];rownames(y)<-x$accession_id
   for(nm in names(y))if(is.character(y[[nm]]))y[[nm]]<-factor(y[[nm]])
   selected<-tryCatch({
    ph<-corehunter::phenotypes(y)
    ch<-corehunter::sampleCore(ph,corehunter::objective('EN','GD'),size=n,
     mode='fast',indices=FALSE,verbose=FALSE)
    ids<-as.character(ch)
    shiny::validate(shiny::need(length(ids)==n&&all(ids%in%x$accession_id),
     'Core Hunter returned unexpected IDs. Check package input conventions.'))
    match(ids,x$accession_id)
   },error=function(e)shiny::validate(shiny::need(FALSE,
    paste('Core Hunter could not run:',conditionMessage(e)))))
   # Core Hunter is unstratified: show observed group counts rather than requested allocation.
   actual<-table(factor(x$group[selected],levels=allocation$group))
   allocation$selected<-as.integer(actual)
  } else if(method=='samplecore'){
   shiny::validate(shiny::need(requireNamespace('SampleCore',quietly=TRUE),
    'Install the supplied SampleCore package and its dependencies first.'))
   # Use exact allocation counts and published SampleCore random selection.
   selected<-tryCatch({
    picked<-SampleCore::select.random(data=x,names='accession_id',group='group',alloc=alloc)
    ids<-unlist(picked,use.names=FALSE)
    shiny::validate(shiny::need(length(ids)==n&&all(ids%in%x$accession_id)&&
     !anyDuplicated(ids),'SampleCore returned invalid or duplicate IDs'))
    match(ids,x$accession_id)
   },error=function(e)shiny::validate(shiny::need(FALSE,
    paste('SampleCore could not run:',conditionMessage(e)))))
  } else {
   ix<-split(seq_len(nrow(x)),x$group)
   selected<-unlist(lapply(names(alloc),function(g){
    k<-alloc[[g]];members<-ix[[g]]
    if(k==0)return(integer(0))
    if(method=='maximin')choose_maximin(members,k,D) else sample(members,k)
   }),use.names=FALSE)
  }
  selected<-as.integer(selected)
  shiny::validate(shiny::need(length(selected)==n&&!anyDuplicated(selected),
   'Selection size or uniqueness validation failed'))
  chosen<-x[selected,,drop=FALSE]
  # Per-descriptor representation metrics, not an EvaluateCore substitute.
  evals<-lapply(traits,function(t){
   a<-x[[t]];b<-chosen[[t]]
   if(is.numeric(a)){
    full_range<-diff(range(a));core_range<-diff(range(b))
    data.frame(descriptor=t,type='numeric',
     full_mean=mean(a),core_mean=mean(b),
     mean_difference=mean(b)-mean(a),
     range_retention_pct=if(full_range>0)100*core_range/full_range else NA_real_,
     category_retention_pct=NA_real_)
   } else {
    data.frame(descriptor=t,type='categorical',full_mean=NA_real_,core_mean=NA_real_,
     mean_difference=NA_real_,range_retention_pct=NA_real_,
     category_retention_pct=100*length(unique(b))/length(unique(a)))
   }
  })
  list(selected=chosen,allocation=allocation,evaluation=do.call(rbind,evals),
   n=n,method=method,strategy=strategy)
 })
 output$sc_alloc<-shiny::renderTable(result()$allocation,striped=TRUE)
 output$sc_selected<-shiny::renderTable(result()$selected,striped=TRUE)
 output$sc_eval<-shiny::renderTable(result()$evaluation,striped=TRUE,digits=3)
 output$sc_notes<-shiny::renderText({
  z<-result()
  paste0('Selected ',z$n,' accessions using ',z$method,'; allocation: ',z$strategy,'.\n',
   'Native stratified random and farthest-first strategies are PlantTraitQuant implementations.\n',
   'SampleCore mode calls SampleCore::select.random using PlantTraitQuant exact group allocations.\n',
   'Core Hunter mode calls corehunter::phenotypes and sampleCore with EN/GD, fast mode. It does NOT enforce group allocation.\n',
   'Distance-weighted allocation uses within-group mean Gower dissimilarity.\n',
   'Evaluation shows mean difference, numeric range retention, and categorical level retention; it does not establish equivalence of full multivariate distributions.\n',
   'A selected subset should be independently validated before breeding or conservation decisions.')
 })
 for(nm in c('selected','alloc','eval'))local({
  key<-nm;field<-switch(key,selected='selected',alloc='allocation',eval='evaluation')
  output[[paste0('sc_',key,'_csv')]]<-shiny::downloadHandler(
   filename=function()paste0('PlantTraitQuant_stratified_',key,'.csv'),
   content=function(file)utils::write.csv(result()[[field]],file,row.names=FALSE,na=''))
 })
}
