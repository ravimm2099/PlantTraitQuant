
# Core collection quantitative quality assessment, numeric and categorical descriptors.
ptq_core_quality_ui <- function() {
 shiny::tabPanel('Core quality evaluation',
  shiny::h3('Quantitative core and mini-core evaluation'),
  shiny::p('Compare a selected accession subset with its full source collection. Include a collection column labelled Full, Core, or Mini-core, or upload a separate selection CSV.'),
  shiny::fluidRow(
   shiny::column(4,shiny::div(class='well',
    shiny::fileInput('cq_data','Full accession descriptor CSV',accept='.csv'),
    shiny::fileInput('cq_selection','Optional selection CSV',accept='.csv'),
    shiny::selectInput('cq_group','Selection grouping column (when no selection CSV)',
     choices=c('Choose after upload'='')),
    shiny::numericInput('cq_reps','Random subset comparisons',199,min=19,max=999),
    shiny::numericInput('cq_seed','Random seed',123,min=1,max=1000000),
    shiny::actionButton('cq_go','Evaluate collection',class='btn-primary'),
    shiny::downloadButton('cq_example','Download example descriptors'),
    shiny::downloadButton('cq_sel_example','Download example selection'))),
   shiny::column(8,
    shiny::h4('Dataset checks'),shiny::verbatimTextOutput('cq_checks'),
    shiny::h4('Core quality indices'),shiny::tableOutput('cq_summary'),
    shiny::downloadButton('cq_summary_csv','Download indices'),
    shiny::h4('Descriptor-level retention'),shiny::tableOutput('cq_descriptors'),
    shiny::downloadButton('cq_desc_csv','Download descriptor metrics'),
    shiny::h4('Interpretation and assumptions'),shiny::verbatimTextOutput('cq_notes'))))
}
ptq_core_quality_server <- function(input,output,session) {
 output$cq_example <- shiny::downloadHandler(
  filename=function()'core_quality_accessions.csv',
  content=function(file)file.copy(file.path('data','core_quality_accessions.csv'),file,overwrite=TRUE))
 output$cq_sel_example <- shiny::downloadHandler(
  filename=function()'core_quality_selection.csv',
  content=function(file)file.copy(file.path('data','core_quality_selection.csv'),file,overwrite=TRUE))
 read_csv <- function(file) {
  shiny::req(file)
  d <- tryCatch(utils::read.csv(file$datapath,check.names=FALSE,
   stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
   error=function(e)shiny::validate(shiny::need(FALSE,conditionMessage(e))))
  names(d) <- trimws(names(d))
  shiny::validate(shiny::need('accession_id'%in%names(d),'Required accession_id column is missing.'),
   shiny::need(!anyDuplicated(names(d)),'Duplicate column names.'))
  d$accession_id <- trimws(as.character(d$accession_id))
  shiny::validate(shiny::need(!anyNA(d$accession_id) &&
   all(nzchar(d$accession_id)) && !anyDuplicated(d$accession_id),
   'Accession IDs must be unique and nonempty.'))
  d
 }
 full <- shiny::reactive(read_csv(input$cq_data))
 shiny::observeEvent(input$cq_data,{
  d <- full()
  choices <- setdiff(names(d),'accession_id')
  shiny::updateSelectInput(session,'cq_group',choices=c('No grouping column'='',choices),
   selected=if('collection'%in%choices)'collection' else '')
 })
 selection <- shiny::reactive({
  d <- full()
  if(!is.null(input$cq_selection)) {
   s <- read_csv(input$cq_selection)
   shiny::validate(shiny::need('collection'%in%names(s),
    'Selection CSV must contain accession_id and collection.'))
   shiny::validate(shiny::need(all(s$accession_id%in%d$accession_id),
    'Selection contains accessions not in full dataset.'))
   group <- rep('Not selected',nrow(d))
   group[match(s$accession_id,d$accession_id)] <- as.character(s$collection)
  } else {
   shiny::req(input$cq_group)
   shiny::validate(shiny::need(nzchar(input$cq_group),
    'Upload a selection CSV or select a collection grouping column.'))
   group <- as.character(d[[input$cq_group]])
  }
  group[is.na(group)] <- 'Not selected'
  # Accept both existing application labels and simple Core/Mini-core labels.
  mini <- grepl('mini',group,ignore.case=TRUE)
  core <- mini | grepl('core',group,ignore.case=TRUE)
  shiny::validate(shiny::need(any(core),'No core accessions found. Use Core only / Mini-core (also core).'),
   shiny::need(sum(core)<nrow(d),'Core must be a proper subset of the full collection.'))
  list(core=which(core),mini=which(mini),group=group)
 })
 output$cq_checks <- shiny::renderText({
  if(is.null(input$cq_data))return('Upload full descriptor CSV and a selection CSV, or use a collection grouping column.')
  d <- full()
  if(is.null(input$cq_selection) && (is.null(input$cq_group)||!nzchar(input$cq_group)))
   return(paste('Full accessions:',nrow(d),'\nSelect a grouping column or upload a selection CSV.'))
  s <- selection()
  paste('Full accessions:',nrow(d),'\nCore:',length(s$core),
   '\nMini-core:',length(s$mini),
   '\nSelection must refer to the same full accession set.')
 })
 result <- shiny::eventReactive(input$cq_go,{
  d <- full();s <- selection()
  excluded <- c('accession_id',input$cq_group,'collection')
  descriptors <- setdiff(names(d),excluded)
  shiny::validate(shiny::need(length(descriptors)>0,'No descriptors to evaluate.'))
  reps <- as.integer(input$cq_reps);seed <- as.integer(input$cq_seed)
  shiny::validate(shiny::need(is.finite(reps)&&reps>=19&&reps<=999,'Random comparisons must be 19 to 999.'),
   shiny::need(is.finite(seed)&&seed>0,'Seed must be positive.'))
  # Interpret numeric-looking columns as numeric; others as categorical.
  values <- lapply(d[descriptors],function(x) {
   v <- suppressWarnings(as.numeric(as.character(x)))
   if(all(is.na(x)|is.finite(v)))v else as.character(x)
  })
  names(values) <- descriptors
  metric <- function(idx) {
   rows <- list()
   for(nm in descriptors) {
    v <- values[[nm]]
    if(is.numeric(v)) {
     vf <- v[is.finite(v)];vs <- v[idx];vs <- vs[is.finite(vs)]
     if(length(vf)<2 || length(vs)<2 || !is.finite(stats::sd(vf)) ||
        stats::sd(vf)<=0)next
     avgdiff <- 100*abs(mean(vs)-mean(vf))/abs(mean(vf))
     if(abs(mean(vf))<1e-12)avgdiff <- NA_real_
     meanret <- if(is.finite(avgdiff))as.numeric(avgdiff<=20) else NA_real_
     varratio <- if(stats::var(vf)>0)stats::var(vs)/stats::var(vf) else NA_real_
     range_ret <- 100*diff(range(vs))/diff(range(vf))
     rows[[length(rows)+1]] <- data.frame(descriptor=nm,type='numeric',
      mean_difference_percent=avgdiff,mean_within_20pct=meanret,
      variance_ratio=varratio,range_retention_percent=range_ret,
      category_retention_percent=NA_real_)
    } else {
     vf <- unique(v[!is.na(v) & nzchar(v)])
     vs <- unique(v[idx][!is.na(v[idx]) & nzchar(v[idx])])
     if(!length(vf))next
     rows[[length(rows)+1]] <- data.frame(descriptor=nm,type='categorical',
      mean_difference_percent=NA_real_,mean_within_20pct=NA_real_,
      variance_ratio=NA_real_,range_retention_percent=NA_real_,
      category_retention_percent=100*length(intersect(vf,vs))/length(vf))
    }
   }
   if(!length(rows))return(NULL)
   do.call(rbind,rows)
  }
  set.seed(seed)
  indices <- list(core=s$core)
  if(length(s$mini))indices$mini_core <- s$mini
  summaries <- list();detail <- list()
  for(label in names(indices)) {
   idx <- indices[[label]]
   m <- metric(idx)
   shiny::validate(shiny::need(!is.null(m),'No evaluable descriptors.'))
   m$collection <- label
   detail[[label]] <- m
   # CR = % numeric traits with mean difference <= 20% (only valid means).
   cr <- if(any(!is.na(m$mean_within_20pct)))100*mean(m$mean_within_20pct,na.rm=TRUE) else NA_real_
   md <- if(any(is.finite(m$mean_difference_percent)))mean(m$mean_difference_percent,na.rm=TRUE) else NA_real_
   vr <- if(any(is.finite(m$variance_ratio)))mean(m$variance_ratio,na.rm=TRUE) else NA_real_
   rr <- if(any(is.finite(m$range_retention_percent)))mean(m$range_retention_percent,na.rm=TRUE) else NA_real_
   catret <- if(any(is.finite(m$category_retention_percent)))mean(m$category_retention_percent,na.rm=TRUE) else NA_real_
   # Random baseline for mean numeric range retention; descriptive, not inferential.
   baseline <- replicate(reps,{
    z <- metric(sample.int(nrow(d),length(idx)))
    if(is.null(z)||!any(is.finite(z$range_retention_percent)))NA_real_
    else mean(z$range_retention_percent,na.rm=TRUE)
   })
   baseline <- baseline[is.finite(baseline)]
   p <- if(length(baseline)&&is.finite(rr)) (1+sum(baseline>=rr))/(1+length(baseline)) else NA_real_
   summaries[[label]] <- data.frame(collection=label,n_selected=length(idx),
    selected_percent=100*length(idx)/nrow(d),
    CR_percent=cr,MD_percent=md,mean_variance_ratio=vr,
    mean_range_retention_percent=rr,mean_category_retention_percent=catret,
    random_mean_range_retention_percent=if(length(baseline))mean(baseline) else NA_real_,
    random_upper_tail_fraction=p)
  }
  list(summary=do.call(rbind,summaries),details=do.call(rbind,detail),
   reps=reps)
 })
 output$cq_summary <- shiny::renderTable(result()$summary,striped=TRUE,digits=3)
 output$cq_descriptors <- shiny::renderTable(result()$details,striped=TRUE,digits=3)
 output$cq_notes <- shiny::renderText({
  z <- result()
  paste0('CR%: percent of eligible numeric descriptors with |core mean - full mean| / |full mean| <= 20%.\n',
   'MD%: average absolute percent mean difference for numeric descriptors with nonzero full means.\n',
   'Variance ratio: mean of selected/full sample variance ratios (not a standardized variance-difference VD%).\n',
   'Range retention: 100 * selected observed range / full observed range for numeric descriptors.\n',
   'Category retention: percentage of observed categorical levels represented.\n',
   'Random comparison: ',z$reps,' equal-size random subsets; upper-tail fraction is descriptive, not a formal selection-adjusted p-value.\n',
   'No molecular allele retention, Mahalanobis distance, multivariate coverage or rarefaction is calculated here.\n',
   'Indices are sensitive to missingness, sample size and trait scaling. Compare with documented sampling designs.')
 })
 output$cq_summary_csv <- shiny::downloadHandler(
  filename=function()'PlantTraitQuant_core_quality_indices.csv',
  content=function(file)utils::write.csv(result()$summary,file,row.names=FALSE,na=''))
 output$cq_desc_csv <- shiny::downloadHandler(
  filename=function()'PlantTraitQuant_core_descriptor_retention.csv',
  content=function(file)utils::write.csv(result()$details,file,row.names=FALSE,na=''))
}
