# PlantTraitQuant v1.3.0: germplasm multivariate and nested core visualization.
# Uses base R and recommended cluster for mixed-trait Gower only when requested.
ptq_v13_read <- function(path) {
 x <- utils::read.csv(path, check.names=FALSE, stringsAsFactors=FALSE, na.strings=c('', 'NA','NaN'), fileEncoding='UTF-8-BOM')
 names(x) <- trimws(names(x))
 shiny::validate(shiny::need('accession_id' %in% names(x), 'CSV requires accession_id'),
  shiny::need(nrow(x)>=6, 'At least six accessions are required'))
 x$accession_id <- trimws(as.character(x$accession_id))
 shiny::validate(shiny::need(!anyNA(x$accession_id) && all(nzchar(x$accession_id)) &&
  !anyDuplicated(x$accession_id), 'Accession IDs must be unique and nonempty'))
 if(!'group' %in% names(x))x$group <- 'All'
 x$group <- as.character(x$group)
 shiny::validate(shiny::need(!anyNA(x$group) && all(nzchar(x$group)), 'Group labels must be nonempty'))
 x
}
ptq_v13_matrix <- function(x,traits) {
 shiny::validate(shiny::need(length(traits)>=2, 'Select at least two numeric traits'))
 shiny::validate(shiny::need(all(vapply(x[,traits,drop=FALSE],is.numeric,logical(1))), 'PCA and Mahalanobis D2 require numeric traits'),
  shiny::need(!anyNA(x[,traits,drop=FALSE]), 'Selected traits must have no missing observations'))
 m<-as.matrix(x[,traits,drop=FALSE]);storage.mode(m)<-'double'
 ok<-apply(m,2,sd)>0
 shiny::validate(shiny::need(all(ok), paste('Remove constant traits:',paste(colnames(m)[!ok],collapse=', '))))
 m
}
ptq_v13_core <- function(x,m,corefrac,minifrac,seed) {
 set.seed(seed); n<-nrow(x)
 # Core: ensure representation from each group where core size permits, then farthest-first
 # on standardized Euclidean trait space. Mini-core: nested farthest-first in the core.
 z<-scale(m); D<-as.matrix(stats::dist(z))
 pick<-function(candidates,k,initial=integer(0)) {
  chosen<-unique(initial); pool<-setdiff(candidates,chosen)
  while(length(chosen)<k && length(pool)) {
   if(!length(chosen)){
    score<-rowMeans(D[pool,candidates,drop=FALSE]); j<-pool[which.max(score)]
   } else {
    score<-apply(D[pool,chosen,drop=FALSE],1,min); j<-pool[which.max(score)]
   }
   chosen<-c(chosen,j);pool<-setdiff(pool,j)
  }
  chosen
 }
 nc<-max(2L,min(n,round(n*corefrac))); nm<-max(1L,min(nc,round(nc*minifrac)))
 groups<-split(seq_len(n),x$group)
 initial<-integer(0)
 if(nc>=length(groups)){
  initial<-vapply(groups,function(ix) ix[which.max(rowMeans(D[ix,,drop=FALSE]))],integer(1))
 }
 core<-pick(seq_len(n),nc,initial)
 # Mini-core selection is restricted to previously selected core accessions.
 mini<-pick(core,nm)
 list(core=core,mini=mini,nc=nc,nm=nm,D=D)
}
ptq_v13_pca <- function(m) {
 fit<-stats::prcomp(m,center=TRUE,scale.=TRUE)
 scores<-fit$x[,1:2,drop=FALSE]
 list(fit=fit,scores=scores,variance=100*fit$sdev^2/sum(fit$sdev^2))
}
ptq_v13_d2 <- function(m) {
 covmat<-stats::cov(m); ev<-eigen(covmat,symmetric=TRUE,only.values=TRUE)$values
 shiny::validate(shiny::need(all(is.finite(ev)) && min(ev)>max(ev)*1e-9,
  'Mahalanobis D2 unavailable: trait covariance is singular or nearly singular. Remove redundant traits or use Euclidean distance.'))
 as.matrix(stats::dist(m %*% chol(solve(covmat))))^2
}
ptq_v13_ui <- function() {
 shiny::tabPanel('Core Visualization Dashboard',
  shiny::h3('Core and mini-core collection visualization'),
  shiny::p('Upload accession-level quantitative descriptors. Generate a nested, illustrative core and mini-core with a reproducible farthest-first algorithm, or upload your own selected accession IDs. PCA and dendrograms are computed on the original collection.'),
  shiny::fluidRow(shiny::column(4,
   shiny::fileInput('cv_data','Accession descriptors (CSV)',accept='.csv'),
   shiny::downloadButton('cv_example','Download example data'),
   shiny::uiOutput('cv_traits_ui'),
   shiny::radioButtons('cv_source','Collection membership',choices=c('Generate nested subsets'='generate','Upload core/mini-core membership'='upload')),
   shiny::conditionalPanel("input.cv_source == 'upload'",shiny::fileInput('cv_membership','Membership CSV: accession_id, collection',accept='.csv')),
   shiny::conditionalPanel("input.cv_source == 'generate'",
    shiny::sliderInput('cv_corefrac','Core size (% of original)',min=5,max=60,value=20,step=1),
    shiny::sliderInput('cv_minifrac','Mini-core size (% of core)',min=5,max=80,value=25,step=1),
    shiny::numericInput('cv_seed','Random seed',123,min=1)),
   shiny::actionButton('cv_run','Build visualization',class='btn-primary'),
   shiny::hr(),shiny::selectInput('cv_trait','Trait for boxplot',choices=NULL),
   shiny::selectInput('cv_fig','Figure to download',choices=c('PCA','Dendrogram','Trait boxplot','Range retention','Cluster coverage','Distance heatmap')),
   shiny::downloadButton('cv_png','Export PNG (300 dpi)'),
   shiny::downloadButton('cv_pdf','Export PDF'),
   shiny::downloadButton('cv_members_csv','Export accession membership'),
   shiny::downloadButton('cv_metrics_csv','Export retention metrics')),
  shiny::column(8,
   shiny::h4('Collection sizes and PCA explained variance'),shiny::verbatimTextOutput('cv_summary'),
   shiny::tabsetPanel(
    shiny::tabPanel('PCA',shiny::plotOutput('cv_pca',height='450px')),
    shiny::tabPanel('Dendrogram',shiny::plotOutput('cv_dendro',height='500px')),
    shiny::tabPanel('Trait distributions',shiny::plotOutput('cv_box',height='450px')),
    shiny::tabPanel('Range retention',shiny::plotOutput('cv_retention',height='450px')),
    shiny::tabPanel('Cluster coverage',shiny::plotOutput('cv_clusters',height='450px')),
    shiny::tabPanel('Distance heatmap',shiny::plotOutput('cv_heatmap',height='500px')),
    shiny::tabPanel('Metrics',shiny::tableOutput('cv_metrics')),
    shiny::tabPanel('Accessions',shiny::tableOutput('cv_members_table'))))),
  shiny::p('Interpretation: retention percentages describe selected quantitative trait ranges, not independent validation of genetic representativeness. Generated subsets use standardized Euclidean farthest-first, not SampleCore or Core Hunter.'))
}
ptq_v13_server <- function(input,output,session) {
 output$cv_example<-shiny::downloadHandler(filename=function()'stratified_core_example.csv',content=function(file)file.copy(file.path('data','stratified_core_example.csv'),file,overwrite=TRUE))
 dat<-shiny::reactive({shiny::req(input$cv_data);ptq_v13_read(input$cv_data$datapath)})
 output$cv_traits_ui<-shiny::renderUI({x<-dat();numeric<-setdiff(names(x)[vapply(x,is.numeric,logical(1))],c('accession_id','group'));shiny::selectInput('cv_traits','Numeric descriptors',choices=numeric,selected=numeric,multiple=TRUE)})
 shiny::observeEvent(input$cv_traits,{shiny::updateSelectInput(session,'cv_trait',choices=input$cv_traits,selected=input$cv_traits[1])})
 res<-shiny::eventReactive(input$cv_run,{
  x<-dat();traits<-input$cv_traits;m<-ptq_v13_matrix(x,traits);p<-ptq_v13_pca(m)
  if(identical(input$cv_source,'upload')){
   shiny::req(input$cv_membership)
   mem<-utils::read.csv(input$cv_membership$datapath,stringsAsFactors=FALSE,check.names=FALSE)
   shiny::validate(shiny::need(all(c('accession_id','collection') %in% names(mem)),'Membership CSV needs accession_id and collection'),
    shiny::need(!anyDuplicated(mem$accession_id),'Membership accession IDs must be unique'),
    shiny::need(all(mem$accession_id %in% x$accession_id),'Membership contains unknown accession IDs'))
   values<-tolower(trimws(as.character(mem$collection)))
   shiny::validate(shiny::need(all(values %in% c('core','mini-core','mini_core','minicore')),'Use collection = core or mini-core'))
   core<-match(mem$accession_id[values %in% c('core','mini-core','mini_core','minicore')],x$accession_id)
   mini<-match(mem$accession_id[values %in% c('mini-core','mini_core','minicore')],x$accession_id)
   shiny::validate(shiny::need(length(core)>=2 && length(mini)>=1,'Upload at least two core accessions and one mini-core accession'))
  }else{
   selection<-ptq_v13_core(x,m,input$cv_corefrac/100,input$cv_minifrac/100,as.integer(input$cv_seed))
   core<-selection$core;mini<-selection$mini
  }
  status<-rep('Original only',nrow(x));status[core]<-'Core';status[mini]<-'Mini-core'
  # For dendrogram and cluster coverage, use the same standardized Euclidean distance.
  d<-stats::dist(scale(m));hc<-stats::hclust(d,method='ward.D2')
  k<-min(5L,max(2L,floor(sqrt(nrow(x)/2))))
  clusters<-stats::cutree(hc,k=k)
  metrics<-do.call(rbind,lapply(seq_len(ncol(m)),function(j){
   full<-m[,j];cvals<-m[core,j];mv<-m[mini,j];r<-diff(range(full))
   data.frame(trait=colnames(m)[j],original_mean=mean(full),core_mean=mean(cvals),mini_mean=mean(mv),
    core_range_retention_pct=if(r>0)100*diff(range(cvals))/r else NA_real_,
    mini_range_retention_pct=if(r>0)100*diff(range(mv))/r else NA_real_)
  }))
  members<-data.frame(accession_id=x$accession_id,group=x$group,collection=status,cluster=clusters,
   PC1=p$scores[,1],PC2=p$scores[,2],stringsAsFactors=FALSE)
  list(x=x,m=m,p=p,core=core,mini=mini,status=status,hc=hc,clusters=clusters,metrics=metrics,members=members,d=as.matrix(d))
 })
 colors<-c('Original only'='#8c98a6','Core'='#187b75','Mini-core'='#d28b24')
 plot_pca<-function(z){
  a<-z$p$scores;cols<-unname(colors[z$status]);pch<-c('Original only'=1,'Core'=16,'Mini-core'=17)[z$status]
  plot(a[,1],a[,2],col=cols,pch=pch,cex=1.25,xlab=sprintf('PC1 (%.1f%%)',z$p$variance[1]),ylab=sprintf('PC2 (%.1f%%)',z$p$variance[2]),main='Germplasm trait PCA: original, core and mini-core')
  abline(h=0,v=0,col='grey85',lty=3);legend('topright',legend=names(colors),col=colors,pch=c(1,16,17),bty='n',cex=.85)
 }
 plot_dendro<-function(z){
  plot(z$hc,labels=FALSE,main='Hierarchical clustering (Ward.D2, standardized traits)',xlab='Accessions',sub='')
  ord<-z$hc$order;usr<-par('usr');y<-usr[3]+.03*diff(usr[3:4]);points(seq_along(ord),rep(y,length(ord)),pch=15,cex=.6,col=colors[z$status[ord]],xpd=NA)
  legend('topright',legend=names(colors),col=colors,pch=15,bty='n',cex=.75)
 }
 plot_box<-function(z){
  t<-input$cv_trait;shiny::validate(shiny::need(!is.null(t)&&t%in%colnames(z$m),'Select a numeric trait'))
  v<-z$m[,t];groups<-list(Original=v,Core=v[z$core],`Mini-core`=v[z$mini]);boxplot(groups,col=c('#cbd5e1','#86c8bc','#f2cb8d'),ylab=t,main=paste('Trait distribution:',t));stripchart(groups,vertical=TRUE,method='jitter',add=TRUE,pch=16,cex=.45,col='#344054')
 }
 plot_retention<-function(z){
  a<-t(as.matrix(z$metrics[,c('core_range_retention_pct','mini_range_retention_pct')]));colnames(a)<-z$metrics$trait
  barplot(a,beside=TRUE,col=c('#187b75','#d28b24'),ylim=c(0,110),las=2,ylab='Range retention (%)',main='Trait range coverage vs original');abline(h=100,lty=2,col='grey50');legend('topright',legend=c('Core','Mini-core'),fill=c('#187b75','#d28b24'),bty='n')
 }
 plot_clusters<-function(z){
  k<-sort(unique(z$clusters));tab<-rbind(Original=tabulate(z$clusters,nbins=max(k))[k],Core=tabulate(z$clusters[z$core],nbins=max(k))[k],`Mini-core`=tabulate(z$clusters[z$mini],nbins=max(k))[k]);colnames(tab)<-paste0('Cluster ',k)
  barplot(tab,beside=TRUE,col=c('#8c98a6','#187b75','#d28b24'),ylab='Number of accessions',main='Coverage of original germplasm clusters');legend('topright',legend=rownames(tab),fill=c('#8c98a6','#187b75','#d28b24'),bty='n')
 }
 plot_heatmap<-function(z){
  ord<-z$hc$order;dd<-z$d[ord,ord,drop=FALSE];graphics::image(seq_len(nrow(dd)),seq_len(ncol(dd)),dd,col=grDevices::colorRampPalette(c('#f7fbff','#6baed6','#08306b'))(50),axes=FALSE,xlab='Accessions ordered by clustering',ylab='Accessions',main='Standardized Euclidean distance heatmap');box()
 }
 figs<-list(PCA=plot_pca,Dendrogram=plot_dendro,`Trait boxplot`=plot_box,`Range retention`=plot_retention,`Cluster coverage`=plot_clusters,`Distance heatmap`=plot_heatmap)
 for(key in names(figs))local({nm<-key;fun<-figs[[key]];out<-switch(nm,PCA='cv_pca',Dendrogram='cv_dendro',`Trait boxplot`='cv_box',`Range retention`='cv_retention',`Cluster coverage`='cv_clusters',`Distance heatmap`='cv_heatmap');output[[out]]<-shiny::renderPlot(fun(res()))})
 output$cv_summary<-shiny::renderText({z<-res();sprintf('Original: %d | Core: %d (%.1f%%) | Mini-core: %d (%.1f%% of core)\nPCA variance: PC1 %.1f%%; PC2 %.1f%%',nrow(z$x),length(z$core),100*length(z$core)/nrow(z$x),length(z$mini),100*length(z$mini)/length(z$core),z$p$variance[1],z$p$variance[2])})
 output$cv_metrics<-shiny::renderTable(res()$metrics,digits=2)
 output$cv_members_table<-shiny::renderTable(utils::head(res()$members,100))
 output$cv_members_csv<-shiny::downloadHandler(filename=function()'PTQ_v13_collection_membership.csv',content=function(file)utils::write.csv(res()$members,file,row.names=FALSE))
 output$cv_metrics_csv<-shiny::downloadHandler(filename=function()'PTQ_v13_trait_retention.csv',content=function(file)utils::write.csv(res()$metrics,file,row.names=FALSE))
 output$cv_png<-shiny::downloadHandler(filename=function()paste0('PTQ_',gsub(' ','_',input$cv_fig),'.png'),content=function(file){grDevices::png(file,width=2400,height=1800,res=300);on.exit(grDevices::dev.off());figs[[input$cv_fig]](res())})
 output$cv_pdf<-shiny::downloadHandler(filename=function()paste0('PTQ_',gsub(' ','_',input$cv_fig),'.pdf'),content=function(file){grDevices::pdf(file,width=10,height=7);on.exit(grDevices::dev.off());figs[[input$cv_fig]](res())})
}
ptq_v13_multivar_ui<-function(){
 shiny::tabPanel('Germplasm Multivariate Analysis',
  shiny::h3('PCA, hierarchical clustering and Mahalanobis D2'),
  shiny::p('Upload an accession descriptor CSV. Analyses use selected quantitative traits; PCA and Euclidean clustering standardize variables. Mahalanobis D2 uses the original covariance matrix and requires full rank.'),
  shiny::fluidRow(shiny::column(4,shiny::fileInput('mv_data','Descriptor CSV',accept='.csv'),shiny::downloadButton('mv_example','Download example'),shiny::uiOutput('mv_traits_ui'),shiny::selectInput('mv_dist','Clustering distance',choices=c('Standardized Euclidean'='euclidean','Mahalanobis D2'='mahalanobis')),shiny::numericInput('mv_k','Number of clusters',4,min=2,max=20),shiny::actionButton('mv_run','Run multivariate analysis',class='btn-primary'),shiny::downloadButton('mv_scores_csv','Download PCA scores and clusters')),
  shiny::column(8,shiny::h4('PCA biplot'),shiny::plotOutput('mv_biplot',height='430px'),shiny::h4('Cluster dendrogram'),shiny::plotOutput('mv_dendro',height='430px'),shiny::h4('Explained variance'),shiny::tableOutput('mv_variance'),shiny::h4('Cluster-wise trait means'),shiny::tableOutput('mv_means'))))
}
ptq_v13_multivar_server<-function(input,output,session){
 output$mv_example<-shiny::downloadHandler(filename=function()'stratified_core_example.csv',content=function(file)file.copy(file.path('data','stratified_core_example.csv'),file,overwrite=TRUE))
 dat<-shiny::reactive({shiny::req(input$mv_data);ptq_v13_read(input$mv_data$datapath)})
 output$mv_traits_ui<-shiny::renderUI({x<-dat();a<-setdiff(names(x)[vapply(x,is.numeric,logical(1))],c('accession_id','group'));shiny::selectInput('mv_traits','Numeric traits',choices=a,selected=a,multiple=TRUE)})
 z<-shiny::eventReactive(input$mv_run,{
  x<-dat();m<-ptq_v13_matrix(x,input$mv_traits);p<-ptq_v13_pca(m)
  if(input$mv_dist=='mahalanobis'){
   d2<-ptq_v13_d2(m);d<-stats::as.dist(sqrt(pmax(d2,0)))
  }else d<-stats::dist(scale(m))
  hc<-stats::hclust(d,method='ward.D2');k<-as.integer(input$mv_k)
  shiny::validate(shiny::need(is.finite(k)&&k>=2&&k<=nrow(m),'Clusters must be between 2 and the number of accessions'))
  cluster<-stats::cutree(hc,k=k)
  means<-aggregate(as.data.frame(m),by=list(cluster=cluster),FUN=mean)
  list(x=x,m=m,p=p,hc=hc,cluster=cluster,means=means)
 })
 output$mv_biplot<-shiny::renderPlot({a<-z();s<-a$p$scores;plot(s[,1],s[,2],pch=19,col=a$cluster+1,xlab=sprintf('PC1 (%.1f%%)',a$p$variance[1]),ylab=sprintf('PC2 (%.1f%%)',a$p$variance[2]),main='PCA scores and trait loadings');abline(h=0,v=0,col='grey85',lty=3);load<-a$p$fit$rotation[,1:2,drop=FALSE];mult<-min(diff(range(s[,1])),diff(range(s[,2])))*.3;arrows(0,0,load[,1]*mult,load[,2]*mult,length=.08,col='#b91c1c');text(load[,1]*mult,load[,2]*mult,labels=rownames(load),pos=3,cex=.75,col='#b91c1c')})
 output$mv_dendro<-shiny::renderPlot({a<-z();plot(a$hc,labels=FALSE,main='Germplasm hierarchical clustering',xlab='Accessions');rect.hclust(a$hc,k=length(unique(a$cluster)),border=2:6)})
 output$mv_variance<-shiny::renderTable({a<-z();data.frame(PC=paste0('PC',seq_along(a$p$variance)),variance_pct=round(a$p$variance,2),cumulative_pct=round(cumsum(a$p$variance),2))})
 output$mv_means<-shiny::renderTable(z()$means,digits=3)
 output$mv_scores_csv<-shiny::downloadHandler(filename=function()'PTQ_v13_PCA_clusters.csv',content=function(file){a<-z();utils::write.csv(data.frame(accession_id=a$x$accession_id,group=a$x$group,cluster=a$cluster,PC1=a$p$scores[,1],PC2=a$p$scores[,2]),file,row.names=FALSE)})
}
