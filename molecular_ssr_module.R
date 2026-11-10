
# Diploid codominant SSR genotype analysis. Long format, one accession-locus per row.
ptq_ssr_ui <- function() {
  shiny::tabPanel('Molecular SSR & marker core',
    shiny::h3('Diploid SSR diversity and marker-based core selection'),
    shiny::p('Upload a long-format, codominant diploid genotype CSV: accession_id, locus, allele1, allele2. Allele labels are categorical, not continuous measurements.'),
    shiny::fluidRow(
      shiny::column(4,shiny::div(class='well',
        shiny::fileInput('ssr_file','Upload diploid SSR genotypes CSV',accept='.csv'),
        shiny::numericInput('ssr_core_pct','Marker core size (%)',20,min=1,max=100),
        shiny::numericInput('ssr_mini_pct','Marker mini-core size (%)',5,min=1,max=100),
        shiny::numericInput('ssr_seed','Tie-breaking random seed',123,min=1,max=1000000),
        shiny::actionButton('ssr_run','Analyse SSR and select core',class='btn-primary'),
        shiny::downloadButton('ssr_example','Download example SSR genotypes'),
        shiny::p('Incomplete genotypes are excluded locus-wise. Selection maximizes observed allele coverage, with distance-based tie-breaking.'))),
      shiny::column(8,
        shiny::h4('Genotype validation'),shiny::verbatimTextOutput('ssr_checks'),
        shiny::h4('Locus-level diversity'),shiny::tableOutput('ssr_loci'),
        shiny::downloadButton('ssr_loci_csv','Download locus statistics CSV'),
        shiny::h4('Core and mini-core selections'),shiny::tableOutput('ssr_selected'),
        shiny::downloadButton('ssr_selection_csv','Download selection CSV'),
        shiny::h4('Allele retention summary'),shiny::tableOutput('ssr_retention'),
        shiny::downloadButton('ssr_retention_csv','Download retention CSV'),
        shiny::h4('Method and limitations'),shiny::verbatimTextOutput('ssr_notes'))))
}
ptq_ssr_server <- function(input,output,session) {
  output$ssr_example <- shiny::downloadHandler(
    filename=function()'diploid_ssr_example.csv',
    content=function(file)file.copy(file.path('data','diploid_ssr_example.csv'),file,overwrite=TRUE))
  raw_data <- shiny::reactive({
    shiny::req(input$ssr_file)
    d <- tryCatch(utils::read.csv(input$ssr_file$datapath,
      stringsAsFactors=FALSE,check.names=FALSE,na.strings=c('','NA','NaN'),
      colClasses='character',fileEncoding='UTF-8-BOM'),
      error=function(e)shiny::validate(shiny::need(FALSE,paste('CSV error:',conditionMessage(e)))))
    names(d) <- trimws(names(d))
    required <- c('accession_id','locus','allele1','allele2')
    shiny::validate(shiny::need(all(required %in% names(d)),
      'Required columns: accession_id, locus, allele1, allele2.'),
      shiny::need(!anyDuplicated(names(d)),'Duplicate CSV column names.'),
      shiny::need(nrow(d)>0,'No genotype rows.'))
    for(nm in required) {
      d[[nm]] <- trimws(d[[nm]])
      d[[nm]][!is.na(d[[nm]]) & !nzchar(d[[nm]])] <- NA_character_
    }
    shiny::validate(shiny::need(!anyNA(d$accession_id) && !anyNA(d$locus),
      'Accession and locus IDs cannot be blank.'),
      shiny::need(!anyDuplicated(paste(d$accession_id,d$locus,sep='\r')),
      'Each accession-locus combination must appear exactly once.'),
      shiny::need(!any(xor(is.na(d$allele1),is.na(d$allele2))),
      'A diploid genotype must have both alleles or both missing.'),
      shiny::need(length(unique(d$accession_id))>=3,'At least three accessions required.'),
      shiny::need(length(unique(d$locus))>=2,'At least two loci required.'),
      shiny::need(length(unique(d$accession_id))<=400,
      'This interactive module supports up to 400 accessions.'))
    d
  })
  output$ssr_checks <- shiny::renderText({
    if(is.null(input$ssr_file))return('Upload a diploid SSR CSV or download the example.')
    d <- raw_data()
    ids <- unique(d$accession_id);loc <- unique(d$locus)
    paste0('Accessions: ',length(ids),' | Loci: ',length(loc),
      '\nGenotype records: ',nrow(d),
      '\nMissing accession-locus records: ',length(ids)*length(loc)-nrow(d),
      '\nGenotypes explicitly missing both alleles: ',sum(is.na(d$allele1)),
      '\nDuplicate accession-locus records: 0',
      '\nOnly fully observed diploid genotypes contribute to locus diversity.')
  })
  calc <- shiny::eventReactive(input$ssr_run,{
    d <- raw_data()
    ids <- unique(d$accession_id)
    loci <- unique(d$locus)
    n <- length(ids)
    core_pct <- input$ssr_core_pct;mini_pct <- input$ssr_mini_pct
    shiny::validate(shiny::need(is.finite(core_pct) && core_pct>0 && core_pct<=100,
      'Core percentage must be 1 to 100.'),
      shiny::need(is.finite(mini_pct) && mini_pct>0 && mini_pct<=core_pct,
      'Mini-core percentage must be positive and not exceed core percentage.'))
    kc <- min(n,max(1L,ceiling(n*core_pct/100)))
    km <- min(kc,max(1L,ceiling(n*mini_pct/100)))
    alleles <- vector('list',length(loci));names(alleles) <- loci
    genotypes <- vector('list',length(loci));names(genotypes) <- loci
    stats_list <- vector('list',length(loci))
    for(j in seq_along(loci)) {
      x <- d[d$locus==loci[j],,drop=FALSE]
      complete <- !is.na(x$allele1) & !is.na(x$allele2)
      x <- x[complete,,drop=FALSE]
      genotypes[[j]] <- x
      if(nrow(x)==0) {
        stats_list[[j]] <- data.frame(locus=loci[j],n_genotyped=0L,
          allele_count=0L,Ho=NA_real_,He=NA_real_,PIC=NA_real_)
        alleles[[j]] <- character()
        next
      }
      aa <- c(x$allele1,x$allele2)
      counts <- table(aa)
      p <- as.numeric(counts)/length(aa)
      ho <- mean(x$allele1!=x$allele2)
      he <- 1-sum(p^2)
      # Codominant multi-allelic PIC (Botstein formula).
      pic <- 1-sum(p^2)-2*sum(outer(p^2,p^2,'*')[upper.tri(outer(p,p))])
      stats_list[[j]] <- data.frame(locus=loci[j],n_genotyped=nrow(x),
        allele_count=length(counts),Ho=ho,He=he,PIC=pic)
      alleles[[j]] <- paste(loci[j],names(counts),sep='::')
    }
    stats <- do.call(rbind,stats_list)
    all_alleles <- unique(unlist(alleles,use.names=FALSE))
    shiny::validate(shiny::need(length(all_alleles)>0,
      'No observed alleles. Provide complete diploid genotypes.'))
    presence <- matrix(FALSE,nrow=n,ncol=length(all_alleles),
      dimnames=list(ids,all_alleles))
    for(j in seq_along(loci)) {
      x <- genotypes[[j]]
      if(nrow(x)==0)next
      for(i in seq_len(nrow(x))) {
        presence[x$accession_id[i],paste(loci[j],x$allele1[i],sep='::')] <- TRUE
        presence[x$accession_id[i],paste(loci[j],x$allele2[i],sep='::')] <- TRUE
      }
    }
    # Genotype allele-sharing distance for tie-breaking: missing pairs ignored.
    dist <- matrix(0,n,n,dimnames=list(ids,ids))
    for(i in seq_len(n))for(k in seq_len(n))if(i<k) {
      shared_loci <- 0L;sumdiff <- 0
      for(j in seq_along(loci)) {
        x <- genotypes[[j]]
        ii <- match(ids[i],x$accession_id);kk <- match(ids[k],x$accession_id)
        if(is.na(ii)||is.na(kk))next
        ai <- c(x$allele1[ii],x$allele2[ii])
        ak <- c(x$allele1[kk],x$allele2[kk])
        # Multiset overlap: 0, 1 or 2 shared allele copies.
        common <- sum(pmin(table(factor(ai,levels=unique(c(ai,ak)))),
                           table(factor(ak,levels=unique(c(ai,ak))))))
        sumdiff <- sumdiff+(1-common/2)
        shared_loci <- shared_loci+1L
      }
      dist[i,k] <- dist[k,i] <- if(shared_loci>0)sumdiff/shared_loci else 1
    }
    seed <- as.integer(input$ssr_seed)
    shiny::validate(shiny::need(!is.na(seed) && seed>0,'Seed must be positive.'))
    set.seed(seed)
    choose <- function(candidates,k,initial=integer()) {
      selected <- initial
      while(length(selected)<k) {
        pool <- setdiff(candidates,selected)
        covered <- if(length(selected))colSums(presence[selected,,drop=FALSE])>0 else rep(FALSE,ncol(presence))
        new_count <- vapply(pool,function(i)sum(presence[i,] & !covered),numeric(1))
        best <- pool[new_count==max(new_count)]
        if(length(best)>1L && length(selected)>0) {
          nearest <- vapply(best,function(i)min(dist[i,selected]),numeric(1))
          best <- best[nearest==max(nearest)]
        }
        selected <- c(selected,sample(best,1))
      }
      selected
    }
    core <- choose(seq_len(n),kc)
    # Re-select mini within core independently for allele coverage; guaranteed nested.
    mini <- choose(core,km)
    summary <- data.frame(collection=c('Core','Mini-core'),
      accessions=c(kc,km),
      percentage=c(100*kc/n,100*km/n),
      alleles_retained=c(sum(colSums(presence[core,,drop=FALSE])>0),
                         sum(colSums(presence[mini,,drop=FALSE])>0)),
      total_observed_alleles=length(all_alleles))
    summary$allele_retention_percent <- 100*summary$alleles_retained/summary$total_observed_alleles
    selection <- data.frame(accession_id=ids,
      collection=ifelse(seq_len(n)%in%mini,'Mini-core (also core)',
        ifelse(seq_len(n)%in%core,'Core only','Not selected')),
      stringsAsFactors=FALSE)
    list(loci=stats,summary=summary,selection=selection,
      missing_pairs=sum(upper.tri(dist) & dist==1))
  })
  output$ssr_loci <- shiny::renderTable(calc()$loci,striped=TRUE,digits=4)
  output$ssr_selected <- shiny::renderTable(
    calc()$selection[calc()$selection$collection!='Not selected',,drop=FALSE],
    striped=TRUE)
  output$ssr_retention <- shiny::renderTable(calc()$summary,striped=TRUE,digits=3)
  output$ssr_loci_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_SSR_locus_statistics.csv',
    content=function(file)utils::write.csv(calc()$loci,file,row.names=FALSE,na=''))
  output$ssr_selection_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_SSR_core_selection.csv',
    content=function(file)utils::write.csv(calc()$selection,file,row.names=FALSE))
  output$ssr_retention_csv <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_SSR_allele_retention.csv',
    content=function(file)utils::write.csv(calc()$summary,file,row.names=FALSE))
  output$ssr_notes <- shiny::renderText({
    paste0('Assumptions: diploid codominant SSRs, unphased allele pairs, categorical allele labels.\n',
      'Ho = observed heterozygote fraction; He = 1 - sum(p^2); PIC uses the Botstein multi-allelic expression.\n',
      'Core selection greedily maximizes the count of distinct observed alleles, with allele-sharing distance tie-breaks.\n',
      'The mini-core is selected only from the core. This is a heuristic, not globally optimal allele maximization.\n',
      'Missing genotype calls and unrecorded accession-locus pairs do not contribute to locus frequencies.\n',
      'He is uncorrected gene diversity. No Hardy-Weinberg, population structure, F-statistics or SNP dosage analysis is claimed.\n',
      'Allele retention is based on the observed dataset; it is not a guarantee of conservation of rare alleles.')
  })
}
