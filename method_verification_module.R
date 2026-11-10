# PlantTraitQuant v3.3.0: self-contained numerical verification examples.
# These are internal reference cases and invariants, not independent software benchmarking.
ptq_validation_cases <- function() {
  tol <- 1e-9
  checks <- list()
  add <- function(name, observed, expected, explanation) {
    checks[[length(checks)+1L]] <<- data.frame(
      test=name,observed=as.character(signif(observed,9)),
      expected=as.character(signif(expected,9)),
      status=if(is.finite(observed) && abs(observed-expected)<tol) 'PASS' else 'FAIL',
      interpretation=explanation,stringsAsFactors=FALSE)
  }
  p <- c(0.5,0.5)
  add('Shannon: two equally abundant species',-sum(p*log(p)),log(2),
      'H = ln(2) using natural logarithms')
  add('Gini-Simpson: two equally abundant species',1-sum(p^2),0.5,
      '1 minus squared relative abundance sum')
  add('Effective species: equal pair',exp(-sum(p*log(p))),2,
      'Exponentiated Shannon diversity')
  add('Pielou evenness: equal pair',(-sum(p*log(p)))/log(2),1,
      'Equal abundances imply maximum evenness')
  # Jaccard presence/absence: {A,B} vs {B,C} => 1 - 1/3.
  a <- c(TRUE,TRUE,FALSE); b <- c(FALSE,TRUE,TRUE)
  add('Jaccard dissimilarity: one shared of three',
      1-sum(a & b)/sum(a | b),2/3,
      'Presence/absence union and intersection')
  # Bray-Curtis: |2-0|+|0-2| divided by total abundance 4.
  x <- c(2,0); y <- c(0,2)
  add('Bray-Curtis: disjoint abundances',sum(abs(x-y))/sum(x+y),1,
      'Completely distinct abundances yield dissimilarity 1')
  # RaoQ and FDis reference cases for two equally weighted species
  # with distance 1 and coordinates at 0 and 1.
  d <- matrix(c(0,1,1,0),2,2)
  add('RaoQ: two species at unit distance',as.numeric(t(p)%*%d%*%p),0.5,
      'Double sum of weighted trait dissimilarities')
  coords <- c(0,1); center <- sum(p*coords)
  add('FDis: equal pair at 0 and 1',sum(p*abs(coords-center)),0.5,
      'Weighted mean distance to trait centroid')
  # Cross-community expectation should be symmetric but can have
  # nonzero diagonal. Pure A vs pure B has cross expectation 1.
  w <- rbind(c(1,0),c(0,1),c(.5,.5))
  cross <- w%*%d%*%t(w)
  add('Cross-community: pure A vs pure B',cross[1,2],1,
      'Expected between-site trait distance')
  add('Cross-community: mixed site diagonal',cross[3,3],0.5,
      'Within-site expected distance need not be zero')
  add('Cross-community: matrix symmetry',max(abs(cross-t(cross))),0,
      'Cross-community expectation is symmetric for symmetric trait distances')
  # 2D convex hull area for triangle with vertices (0,0), (1,0), (0,1).
  xy <- matrix(c(0,0,1,0,0,1),ncol=2,byrow=TRUE)
  hull <- grDevices::chull(xy[,1],xy[,2]); poly <- xy[c(hull,hull[1]),,drop=FALSE]
  area <- abs(sum(poly[-nrow(poly),1]*poly[-1,2] -
                  poly[-1,1]*poly[-nrow(poly),2]))/2
  add('2D convex hull: right triangle',area,0.5,
      'Area of unit right triangle equals one half')
  do.call(rbind,checks)
}
ptq_validation_ui <- function() {
  shiny::tabPanel('Method verification',
    shiny::h3('Reproducible numerical checks'),
    shiny::p('Run small analytical reference cases to verify key mathematical definitions used throughout PlantTraitQuant. This is a diagnostic aid, not a substitute for independent package-to-package validation.'),
    shiny::fluidRow(
      shiny::column(4,shiny::wellPanel(
        shiny::actionButton('verify_run','Run numerical checks',class='btn-primary'),
        shiny::downloadButton('verify_download','Download verification CSV'),
        shiny::p('The test cases are deterministic, require no uploaded data, and use only base R and shiny.'))),
      shiny::column(8,
        shiny::h4('Verification summary'),shiny::verbatimTextOutput('verify_summary'),
        shiny::h4('Individual checks'),shiny::tableOutput('verify_table'),
        shiny::h4('Scope and limitations'),
        shiny::p('These tests verify known toy-case answers and invariants. They do not verify the entire GUI, all ecological edge cases, FEve/FDiv algorithms, or agreement with external software.'))))
}
ptq_validation_server <- function(input,output,session) {
  results <- shiny::eventReactive(input$verify_run,ptq_validation_cases())
  output$verify_summary <- shiny::renderText({
    z <- results()
    paste('Passed:',sum(z$status=='PASS'),'/',nrow(z),
      '\nFailed:',sum(z$status!='PASS'),
      '\nStatus:',if(all(z$status=='PASS'))'All included checks passed' else 'Review failed checks')
  })
  output$verify_table <- shiny::renderTable(results(),striped=TRUE,spacing='xs')
  output$verify_download <- shiny::downloadHandler(
    filename=function()'PlantTraitQuant_numerical_verification.csv',
    content=function(file)utils::write.csv(results(),file,row.names=FALSE))
}
