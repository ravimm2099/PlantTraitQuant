# PlantTraitQuant 2.2.0 - independent trait-space analysis
# Required packages: shiny (only). Run: shiny::runApp('.')
library(shiny)
source("mixed_trait_module.R")
source("functional_diversity_module.R")
source("fric2d_module.R")
source("feve_module.R")
source("fdiv2d_module.R")
source("dashboard_module.R")
source("beta_diversity_module.R")
source("functional_beta_module.R")
source("pairwise_functional_module.R")
source("method_verification_module.R")
source("independent_benchmark_module.R")
source("trait_tree_beta_module.R")
source("core_collection_module.R")
source("advanced_core_module.R")
source("molecular_ssr_module.R")
source("molecular_snp_module.R")
source("trait_environment_module.R")
source("ade4_rlq_fourthcorner_module.R")
source("core_quality_module.R")
source("functional_structure_module.R")
source("multidimensional_fric_module.R")
source("stratified_core_module.R")
source("agronomic_genetic_module.R")
source("core_visualization_v13_module.R")

ui <- fluidPage(
 tags$head(
  tags$meta(name='viewport',content='width=device-width, initial-scale=1'),
  tags$title('PlantTraitQuant v1.3.0 | Dr. Manish Mathur'),
  tags$style(HTML("
body {background:#f3f6f3;color:#18382d;font-family:'Segoe UI',Arial,sans-serif;margin:0;}
.container-fluid {padding:0;}
.ptq-topbar {background:#073c2c;color:#fff;padding:12px 25px;display:flex;align-items:center;justify-content:space-between;gap:16px;flex-wrap:wrap;}
.ptq-brand {font-size:25px;font-weight:750;letter-spacing:-.5px;color:white;}
.ptq-brand small {font-size:13px;font-weight:400;margin-left:7px;color:#d8eee2;}
.ptq-topnav {display:flex;gap:10px;align-items:center;flex-wrap:wrap;}
.ptq-topnav a {color:white;padding:8px 11px;text-decoration:none;font-weight:600;}
.ptq-topnav a:hover {background:#2d8656;color:white;border-radius:8px;}
.ptq-layout {display:flex;min-height:calc(100vh - 160px);}
.ptq-sidebar {width:240px;flex:0 0 240px;background:linear-gradient(175deg,#063e2d,#012b20);color:#fff;padding:15px 0;max-height:calc(100vh - 75px);overflow-y:auto;position:sticky;top:0;}
.ptq-sidebar h3 {color:#fff;font-size:19px;padding:2px 20px 12px;margin:0;font-weight:700;}
.ptq-side-link {display:block;color:#f5fff8;padding:10px 17px;font-size:13px;text-decoration:none;line-height:1.3;border-left:4px solid transparent;}
.ptq-side-link:hover,.ptq-side-link:focus {background:#34834e;color:#fff;text-decoration:none;border-left-color:#b5e8b0;}
.ptq-nav-icon {display:inline-block;width:23px;font-size:17px;}
.ptq-main {min-width:0;flex:1;padding:18px 20px 28px;}
.ptq-main .nav-tabs {display:none!important;}
.ptq-main .tab-content {background:transparent;}
.ptq-hero {background:linear-gradient(105deg,rgba(255,255,255,.95),rgba(243,248,237,.94)),linear-gradient(150deg,#e7ede5,#fbf1df);padding:34px;border-radius:16px;border:1px solid #dfe9d9;box-shadow:0 6px 20px #163d2511;}
.ptq-hero-grid {display:grid;grid-template-columns:minmax(0,1.7fr) minmax(250px,.85fr);gap:28px;align-items:start;}
.ptq-app-title {font-family:Georgia,serif;font-weight:800;color:#094a25;font-size:clamp(32px,4vw,56px);margin:0 0 9px;}
.ptq-pill {display:inline-block;background:#38752a;color:#fff;border-radius:14px;font:700 18px Arial;padding:7px 13px;vertical-align:middle;}
.ptq-subtitle {font-size:23px;font-weight:600;font-style:italic;color:#5c3824;line-height:1.3;margin-bottom:18px;}
.ptq-description {font-size:16px;color:#1c2f27;max-width:800px;line-height:1.6;}
.ptq-chips {display:flex;flex-wrap:wrap;gap:9px;margin-top:22px;}
.ptq-chip {background:#e7ebd8;border:1px solid #d8e0c7;border-radius:25px;padding:8px 12px;font-size:12px;font-weight:700;color:#27412e;}
.ptq-profile {padding:2px 0;}
.ptq-profile img {width:100%;max-width:245px;height:230px;object-fit:cover;object-position:center 24%;border:5px solid white;border-radius:19px;box-shadow:0 3px 15px #213b2c38;display:block;margin:0 auto 15px;}
.ptq-profile h2 {font-size:23px;color:#083d2a;margin:0 0 4px;border-left:6px solid #267c39;padding-left:10px;}
.ptq-credentials {font-size:13px;font-weight:650;font-style:italic;margin:0 0 9px;}
.ptq-profile p {font-size:13px;line-height:1.4;margin-bottom:7px;}
.ptq-profile a {overflow-wrap:anywhere;}
.ptq-stats {display:grid;grid-template-columns:repeat(4,minmax(0,1fr));background:#fffaf0;border:1px solid #e6dccc;border-radius:13px;margin:18px 0;box-shadow:0 2px 10px #153a2810;}
.ptq-stat {padding:14px 20px;text-align:center;border-right:1px solid #b8c6b3;}
.ptq-stat:last-child {border-right:0;}
.ptq-stat b {font-size:29px;color:#19552b;display:block;line-height:1.1;}
.ptq-stat span {font-size:12px;color:#2b332e;}
.ptq-tile-grid {display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;}
.ptq-tile {display:block;min-height:105px;padding:17px 14px;border-radius:13px;text-decoration:none!important;color:#142f28!important;border:1px solid #ffffffaa;transition:transform .12s;}
.ptq-tile:hover {transform:translateY(-2px);box-shadow:0 4px 14px #2b4e281c;}
.ptq-tile strong {display:block;font-size:16px;margin-bottom:9px;}
.ptq-tile small {font-size:12px;line-height:1.5;}
.t-green{background:#dfeecf}.t-blue{background:#dceaf9}.t-purple{background:#e9dcf6}.t-orange{background:#fae0d3}
.t-rose{background:#f8dfdb}.t-cyan{background:#d8f1ef}.t-gold{background:#fff0c8}.t-lavender{background:#e9e5f4}
.ptq-footer {background:#fff8e8;padding:18px 25px;border-top:2px solid #e5d7bb;display:flex;gap:18px;align-items:center;justify-content:space-between;flex-wrap:wrap;color:#14392b;font-size:13px;}
.ptq-footer b{color:#16552f}.ptq-footer img{height:45px;vertical-align:middle;margin-right:8px;}
.well {background:#fff;border:1px solid #dcebe3;border-radius:10px;}
.btn-primary {background:#176d51;border-color:#176d51;}
.ptq-note {padding:12px;background:#eaf6f0;border-left:4px solid #21875e;margin:10px 0;}
@media(max-width:1150px){.ptq-tile-grid{grid-template-columns:repeat(2,minmax(0,1fr));}.ptq-hero-grid{grid-template-columns:1fr;}.ptq-profile{max-width:600px;}}
@media(max-width:750px){.ptq-layout{display:block}.ptq-sidebar{position:static;width:100%;max-height:200px;overflow:auto}.ptq-main{padding:12px}.ptq-stats{grid-template-columns:repeat(2,1fr)}.ptq-tile-grid{grid-template-columns:1fr}.ptq-hero{padding:20px}.ptq-topbar{padding:10px}}
"))
 ),
 div(class='ptq-topbar',
  div(class='ptq-brand','❧ PlantTraitQuant',tags$small('Plant Trait, Functional Diversity & Germplasm Analysis')),
  div(class='ptq-topnav',tags$a(href='#',onclick="Shiny.setInputValue('ptq_nav', 'Home & guide', {priority:'event'}); return false;",'⌂ Home'),tags$a(href='#',onclick="Shiny.setInputValue('ptq_nav', 'About', {priority:'event'}); return false;",'ⓘ About'),tags$a(href='#',onclick="Shiny.setInputValue('ptq_nav', 'Home & guide', {priority:'event'}); return false;",'▣ User Guide'),tags$a(href='#',onclick="Shiny.setInputValue('ptq_nav', 'Species traits', {priority:'event'}); return false;",'▤ Data'))),
 div(class='ptq-layout',
  div(class='ptq-sidebar',tags$h3('Main Menu'),
   tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Home & guide', {priority:'event'}); return false;", span(class='ptq-nav-icon','⌂'), 'Home'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Interactive Workflow Guide', {priority:'event'}); return false;", span(class='ptq-nav-icon','?'), 'Interactive Workflow Guide'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Species traits', {priority:'event'}); return false;", span(class='ptq-nav-icon','▦'), 'Trait & Data Utilities'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Species abundance', {priority:'event'}); return false;", span(class='ptq-nav-icon','▥'), 'Species Abundance'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Trait space (PCA)', {priority:'event'}); return false;", span(class='ptq-nav-icon','◈'), 'Trait PCA'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Mixed traits (Gower + PCoA)', {priority:'event'}); return false;", span(class='ptq-nav-icon','◇'), 'Mixed Traits & Gower'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Community beta diversity', {priority:'event'}); return false;", span(class='ptq-nav-icon','❀'), 'Taxonomic Beta Diversity'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Community functional diversity', {priority:'event'}); return false;", span(class='ptq-nav-icon','◉'), 'Community Functional Diversity'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Functional richness (2D FRic)', {priority:'event'}); return false;", span(class='ptq-nav-icon','⬡'), 'Functional Richness'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Functional evenness (FEve)', {priority:'event'}); return false;", span(class='ptq-nav-icon','◎'), 'Functional Evenness'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Functional divergence (2D FDiv)', {priority:'event'}); return false;", span(class='ptq-nav-icon','◈'), 'Functional Divergence'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'FD dashboard (all indices)', {priority:'event'}); return false;", span(class='ptq-nav-icon','▤'), 'FD Dashboard'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Functional beta diversity', {priority:'event'}); return false;", span(class='ptq-nav-icon','◇'), 'Functional Beta Diversity'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Pairwise functional differences', {priority:'event'}); return false;", span(class='ptq-nav-icon','↔'), 'Pairwise Functional Differences'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Trait-tree beta partition', {priority:'event'}); return false;", span(class='ptq-nav-icon','♧'), 'Trait-Tree Beta'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Trait-environment associations', {priority:'event'}); return false;", span(class='ptq-nav-icon','⌁'), 'Trait–Environment'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'RLQ & Fourth-Corner', {priority:'event'}); return false;", span(class='ptq-nav-icon','⟐'), 'RLQ & Fourth-Corner'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Functional structure & scales', {priority:'event'}); return false;", span(class='ptq-nav-icon','▧'), 'Functional Structure & Scales'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Multidimensional FRic', {priority:'event'}); return false;", span(class='ptq-nav-icon','⬡'), 'Multidimensional FRic'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Core and mini-core collection', {priority:'event'}); return false;", span(class='ptq-nav-icon','❁'), 'Core & Mini-Core'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Advanced core selection', {priority:'event'}); return false;", span(class='ptq-nav-icon','✥'), 'Advanced Core Selection'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Core quality evaluation', {priority:'event'}); return false;", span(class='ptq-nav-icon','▣'), 'Core Quality'),
   tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Stratified Core Selection', {priority:'event'}); return false;", span(class='ptq-nav-icon','✣'), 'Stratified Core Selection'),
   tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Agronomic & Genetic Variability', {priority:'event'}); return false;", span(class='ptq-nav-icon','▦'), 'Agronomic & Genetic Variability'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Core Visualization Dashboard', {priority:'event'}); return false;", span(class='ptq-nav-icon','◉'), 'Core Visualization Dashboard'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Germplasm Multivariate Analysis', {priority:'event'}); return false;", span(class='ptq-nav-icon','◇'), 'Germplasm Multivariate'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Molecular SSR & marker core', {priority:'event'}); return false;", span(class='ptq-nav-icon','⌘'), 'SSR Molecular Diversity'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'SNP diversity & genomic core', {priority:'event'}); return false;", span(class='ptq-nav-icon','⌘'), 'SNP Molecular Diversity'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Method verification', {priority:'event'}); return false;", span(class='ptq-nav-icon','✓'), 'Method Verification'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Independent benchmarking', {priority:'event'}); return false;", span(class='ptq-nav-icon','▥'), 'Independent Benchmarking'),
tags$a(class='ptq-side-link', href='#', onclick="Shiny.setInputValue('ptq_nav', 'About', {priority:'event'}); return false;", span(class='ptq-nav-icon','ⓘ'), 'About')),
  div(class='ptq-main',
   tabsetPanel(id='tabs',type='tabs',
    tabPanel('Home & guide',
     div(class='ptq-hero',
      div(class='ptq-hero-grid',
       div(
        div(class='ptq-app-title','❧ PlantTraitQuant'),
        div(class='ptq-subtitle','An Integrated R Shiny Toolkit for Plant Trait, Functional Diversity and Germplasm Quantitative Analysis'),
        p(class='ptq-description','An open-source R Shiny application for quantitative analysis of plant traits, taxonomic and functional diversity, trait–environment relationships, community structure, core collection evaluation and molecular diversity, with downloadable results.'),
        div(class='ptq-chips',
         span(class='ptq-chip','▥ Trait Analysis'),span(class='ptq-chip','❀ Functional Diversity'),
         span(class='ptq-chip','✥ Germplasm Analysis'),span(class='ptq-chip','⌖ Spatial & Multiscale'))
       ),
       div(class='ptq-profile',
        tags$img(src='developer_portrait.jpg',alt='Portrait supplied by Dr. Manish Mathur'),
        tags$h2('Dr. Manish Mathur'),
        p(class='ptq-credentials','Ph.D. (Plant Ecology), ARS-NET, Crop Protection'),
        p('ICAR-Central Arid Zone Research Institute (ICAR-CAZRI), Jodhpur, Rajasthan, India'),
        p(tags$a(href='mailto:manish.mathur@icar.org.in','manish.mathur@icar.org.in')),
        p(tags$a(href='mailto:eco5320@gmail.com','eco5320@gmail.com'))
       )
      )
     ),
     div(class='ptq-stats',
      div(class='ptq-stat',tags$b('Multiple'),'Quantitative methods'),
      div(class='ptq-stat',tags$b('24'),'Analytical module files'),
      div(class='ptq-stat',tags$b('2–5D'),'Multidimensional FRic'),
      div(class='ptq-stat',tags$b('CSV'),'Upload & download')),
     div(class='ptq-tile-grid',
      tags$a(class='ptq-tile t-green', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Community beta diversity', {priority:'event'}); return false;", tags$strong('Taxonomic Diversity'),tags$small('Alpha, beta and community indices')),
tags$a(class='ptq-tile t-blue', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Community functional diversity', {priority:'event'}); return false;", tags$strong('Functional Diversity'),tags$small('FDis, RaoQ, FRic, FEve, FDiv')),
tags$a(class='ptq-tile t-purple', href='#', onclick="Shiny.setInputValue('ptq_nav', 'RLQ & Fourth-Corner', {priority:'event'}); return false;", tags$strong('Trait–Environment'),tags$small('RLQ and fourth-corner analysis')),
tags$a(class='ptq-tile t-orange', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Functional structure & scales', {priority:'event'}); return false;", tags$strong('Functional Structure'),tags$small('Originality, redundancy and specialization')),
tags$a(class='ptq-tile t-rose', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Multidimensional FRic', {priority:'event'}); return false;", tags$strong('Multidimensional FRic'),tags$small('Convex-hull richness in 2–5 dimensions')),
tags$a(class='ptq-tile t-cyan', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Core quality evaluation', {priority:'event'}); return false;", tags$strong('Core & Mini-Core'),tags$small('Accession selection and quality evaluation')),
tags$a(class='ptq-tile t-gold', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Molecular SSR & marker core', {priority:'event'}); return false;", tags$strong('Molecular Diversity'),tags$small('SSR and SNP marker analyses')),
tags$a(class='ptq-tile t-lavender', href='#', onclick="Shiny.setInputValue('ptq_nav', 'Functional structure & scales', {priority:'event'}); return false;", tags$strong('Spatial & Multiscale'),tags$small('SES, null models and site accumulation'))
     ),
     tags$p(style='font-size:12px;color:#53695b;margin-top:16px',
      'Scientific note: some measures are exploratory proxies; multidimensional convex-hull FRic is not a probabilistic hypervolume. See INDEX_REGISTRY.md.')
    ),
    tabPanel('Interactive Workflow Guide', tags$iframe(src='interactive_workflow.html',style='width:100%;height:84vh;border:0;',title='Interactive workflow guide')),
    tabPanel('Species traits',
      fluidRow(
        column(4,div(class='well',
          fileInput('trait_file','Upload species trait CSV',accept=c('.csv','text/csv')),
          selectInput('missing_rule','Missing values',choices=c('Use available values for each trait'='available','Complete observations only'='complete')),
          checkboxGroupInput('calculations','Choose calculations',choices=c('Trait summary'='summary','Species summary'='species','Trait correlations'='correlation'),selected='summary'),
          actionButton('run','Calculate',class='btn-primary'),
          downloadButton('example','Download example CSV'))),
        column(8,
          h4('Validation'),verbatimTextOutput('validation'),
          h4('Data preview'),tableOutput('preview'),
          h4('Results'),tableOutput('results'),
          downloadButton('download','Download calculated results CSV')))),
    tabPanel('Species abundance',
      h3('Independent species abundance analysis'),
      p('Upload a site-by-species CSV. The first column must be site_id; remaining columns contain nonnegative counts, cover, or other comparable abundance measures. No trait file is needed.'),
      fluidRow(
        column(4, div(class='well',
          fileInput('abundance_file','Upload site abundance CSV',accept=c('.csv','text/csv')),
          selectInput('abundance_missing','Missing abundance cells',
            choices=c('Stop and report missing values'='stop','Treat missing as zero (only if justified)'='zero')),
          actionButton('run_abundance','Calculate diversity',class='btn-primary'),
          downloadButton('abundance_example','Download example abundance CSV'))),
        column(8,
          h4('Validation and data quality'),verbatimTextOutput('abundance_validation'),
          h4('Data preview'),tableOutput('abundance_preview'),
          h4('Diversity results'),tableOutput('abundance_results'),
          downloadButton('abundance_download','Download abundance results CSV')))),
    tabPanel('Trait space (PCA)',
      h3('Species trait-space analysis'),
      p('Upload a species trait CSV independently. Choose numeric traits for standardized PCA. No abundance file or extra R package is required.'),
      fluidRow(
        column(4,div(class='well',
          fileInput('pca_file','Upload species trait CSV',accept=c('.csv','text/csv')),
          uiOutput('pca_trait_choices'),
          actionButton('run_pca','Calculate PCA',class='btn-primary'),
          downloadButton('pca_example','Download example trait CSV'),
          p('Rows with missing selected traits are excluded only from PCA. Constant traits are reported and excluded.'))),
        column(8,
          h4('Data validation'),verbatimTextOutput('pca_validation'),
          h4('Species ordination (PC1 versus PC2)'),plotOutput('pca_plot',height='480px'),
          h4('Explained variance'),tableOutput('pca_variance'),
          h4('Species scores'),tableOutput('pca_scores'),
          downloadButton('pca_scores_download','Download PCA scores CSV'),
          h4('Trait loadings'),tableOutput('pca_loadings'),
          downloadButton('pca_loadings_download','Download trait loadings CSV')))),
    ptq_mixed_ui(),
    ptq_fd_ui(),
    ptq_fric_ui(),
    ptq_feve_ui(),
    ptq_fdiv_ui(),
    ptq_dashboard_ui(),
    ptq_beta_ui(),
    ptq_fbeta_ui(),
    ptq_pairwise_ui(),
    ptq_validation_ui(),
    ptq_benchmark_ui(),
    ptq_treebeta_ui(),
    ptq_core_ui(),
    ptq_advanced_core_ui(),
    ptq_ssr_ui(),
    ptq_snp_ui(),
    ptq_te_ui(),
    ptq_ade4_ui(),
    ptq_core_quality_ui(),
    ptq_structure_ui(),
    ptq_hull_ui(),
    ptq_stratified_ui(),
    ptq_agronomic_ui(),
    ptq_v13_ui(),
    ptq_v13_multivar_ui(),
    tabPanel('About',
     h2('About PlantTraitQuant'),
     div(class='ptq-profile',tags$img(src='developer_portrait.jpg',alt='Developer portrait'),
      tags$h2('Dr. Manish Mathur'),
      p(class='ptq-credentials','Ph.D. (Plant Ecology), ARS-NET, Crop Protection'),
      p('ICAR-Central Arid Zone Research Institute (ICAR-CAZRI), Jodhpur, Rajasthan, India'),
      p(tags$a(href='mailto:manish.mathur@icar.org.in','manish.mathur@icar.org.in')),
      p(tags$a(href='mailto:eco5320@gmail.com','eco5320@gmail.com'))),
     h3('Scientific scope and limitations'),
     p('The application integrates quantitative methods for plant traits, community and functional diversity, molecular markers, and germplasm core selection. Some methods are exploratory proxies and not independently benchmarked.'),
     p('Multidimensional convex-hull FRic is not a probabilistic kernel hypervolume. Consult INDEX_REGISTRY.md and README.md for definitions and limitations.'),
     h3('Citation'),p('Mathur, M. (2026). PlantTraitQuant (software release) [Computer software]. See CITATION.cff.'))
   )
  )
 ),
 div(class='ptq-footer',
  div(tags$b('ICAR-Central Arid Zone Research Institute (ICAR-CAZRI)'),tags$div('Jodhpur, Rajasthan, India')),
  div('Developed by: ',tags$b('Dr. Manish Mathur'),tags$div('© 2026 PlantTraitQuant | R Shiny')))
)

server <- function(input,output,session) {
  observeEvent(input$ptq_nav, {
    req(input$ptq_nav)
    updateTabsetPanel(session,'tabs',selected=input$ptq_nav)
  })
  ptq_mixed_server(input,output,session)
  ptq_fd_server(input,output,session)
  ptq_fric_server(input,output,session)
  ptq_feve_server(input,output,session)
  ptq_fdiv_server(input,output,session)
  ptq_dashboard_server(input,output,session)
  ptq_beta_server(input,output,session)
  ptq_fbeta_server(input,output,session)
  ptq_pairwise_server(input,output,session)
  ptq_validation_server(input,output,session)
  ptq_benchmark_server(input,output,session)
  ptq_treebeta_server(input,output,session)
  ptq_core_server(input,output,session)
  ptq_advanced_core_server(input,output,session)
  ptq_ssr_server(input,output,session)
  ptq_snp_server(input,output,session)
  ptq_te_server(input,output,session)
  ptq_ade4_server(input,output,session)
  ptq_core_quality_server(input,output,session)
  ptq_structure_server(input,output,session)
  ptq_hull_server(input,output,session)
  ptq_stratified_server(input,output,session)
  ptq_agronomic_server(input,output,session)
  ptq_v13_server(input,output,session)
  ptq_v13_multivar_server(input,output,session)
  example_path <- file.path('data','species_traits_example.csv')
  output$example <- downloadHandler(filename=function() 'species_traits_example.csv',content=function(file) file.copy(example_path,file))
  raw_data <- reactive({
    req(input$trait_file)
    tryCatch(read.csv(input$trait_file$datapath,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM'),
             error=function(e) validate(need(FALSE,paste('CSV could not be read:',conditionMessage(e)))))
  })
  cleaned <- reactive({
    d <- raw_data()
    names(d) <- trimws(names(d))
    validate(need('species_id' %in% names(d),'Missing required column: species_id. Download the example CSV for the correct format.'))
    d$species_id <- trimws(as.character(d$species_id))
    validate(need(nrow(d)>0,'The uploaded file has no rows.'),
             need(all(nzchar(d$species_id)),'Some species_id entries are blank.'),
             need(!anyDuplicated(d$species_id),'Duplicate species_id values found. Use unique species identifiers.'))
    traits <- setdiff(names(d),'species_id')
    for (nm in traits) {
      if (is.character(d[[nm]])) {
        x <- suppressWarnings(as.numeric(d[[nm]]))
        if (all(is.na(d[[nm]]) | !is.na(x))) d[[nm]] <- x
      }
    }
    numeric_traits <- traits[vapply(d[traits],is.numeric,logical(1))]
    validate(need(length(numeric_traits)>0,'At least one numeric trait column is required.'))
    list(data=d,traits=numeric_traits)
  })
  output$validation <- renderText({
    if (is.null(input$trait_file)) return('No CSV uploaded. Download the example file or choose your own CSV.')
    x <- cleaned(); d <- x$data
    paste0('Rows: ',nrow(d),' | Columns: ',ncol(d),' | Numeric traits: ',length(x$traits),
           '\nMissing numeric values: ',sum(is.na(d[x$traits])),
           '\nNon-numeric columns (excluded from numerical calculations): ',
           paste(setdiff(names(d),c('species_id',x$traits)),collapse=', '),
           '\nIdentifiers: unique; data ready for selected analyses.')
  })
  output$preview <- renderTable({head(cleaned()$data,10)},striped=TRUE,spacing='xs')
  computed <- eventReactive(input$run,{
    x <- cleaned(); d <- x$data; tr <- x$traits
    validate(need(length(input$calculations)>0,'Select at least one calculation.'))
    if (identical(input$missing_rule,'complete')) d <- d[complete.cases(d[tr]),,drop=FALSE]
    validate(need(nrow(d)>0,'No observations remain after the missing-value rule.'))
    parts <- list()
    if ('summary' %in% input$calculations) {
      parts[['summary']] <- do.call(rbind,lapply(tr,function(nm) {
        v <- d[[nm]];v <- v[is.finite(v)]
        data.frame(analysis='trait_summary',variable=nm,species_id='',n=length(v),
                   mean=if(length(v))mean(v) else NA_real_,
                   sd=if(length(v)>1)sd(v) else NA_real_,
                   min=if(length(v))min(v) else NA_real_,
                   max=if(length(v))max(v) else NA_real_,
                   value=NA_real_)
      }))
    }
    if ('species' %in% input$calculations) {
      parts[['species']] <- do.call(rbind,lapply(seq_len(nrow(d)),function(i) {
        v <- unlist(d[i,tr,drop=FALSE],use.names=FALSE);v <- v[is.finite(v)]
        data.frame(analysis='species_summary',variable='all_numeric_traits',species_id=d$species_id[i],n=length(v),
                   mean=if(length(v))mean(v) else NA_real_,
                   sd=if(length(v)>1)sd(v) else NA_real_,
                   min=if(length(v))min(v) else NA_real_,
                   max=if(length(v))max(v) else NA_real_,value=NA_real_)
      }))
    }
    if ('correlation' %in% input$calculations) {
      validate(need(length(tr)>=2,'Correlations require at least two numeric traits.'))
      m <- suppressWarnings(cor(d[tr],use='pairwise.complete.obs'))
      parts[['correlation']] <- do.call(rbind,lapply(seq_len(nrow(m)),function(i)
        data.frame(analysis='trait_correlation',variable=rownames(m)[i],species_id=colnames(m),
                   n=NA_integer_,mean=NA_real_,sd=NA_real_,min=NA_real_,max=NA_real_,value=as.numeric(m[i,]))))
    }
    do.call(rbind,parts)
  })
  output$results <- renderTable({head(computed(),100)},striped=TRUE,spacing='xs')
  abundance_example_path <- file.path('data','site_abundance_example.csv')
  output$abundance_example <- downloadHandler(
    filename=function() 'site_abundance_example.csv',
    content=function(file) file.copy(abundance_example_path,file))
  abundance_data <- reactive({
    req(input$abundance_file)
    d <- tryCatch(read.csv(input$abundance_file$datapath, check.names=FALSE,
                  stringsAsFactors=FALSE, na.strings=c('','NA','NaN'),
                  fileEncoding='UTF-8-BOM'),
                  error=function(e) validate(need(FALSE,paste('CSV could not be read:',conditionMessage(e)))))
    names(d) <- trimws(names(d))
    validate(need(nrow(d)>0,'The file contains no sites.'),
             need('site_id' %in% names(d),'Required column site_id is missing.'),
             need(!anyDuplicated(names(d)),'Duplicate column names found.'))
    d$site_id <- trimws(as.character(d$site_id))
    validate(need(all(!is.na(d$site_id) & nzchar(d$site_id)),'Blank site_id found.'),
             need(!anyDuplicated(d$site_id),'site_id values must be unique.'))
    species <- setdiff(names(d),'site_id')
    validate(need(length(species)>0,'Add at least one species abundance column.'))
    for (nm in species) {
      values <- suppressWarnings(as.numeric(d[[nm]]))
      invalid <- !is.na(d[[nm]]) & is.na(values)
      validate(need(!any(invalid),paste('Non-numeric abundance in:',nm)))
      d[[nm]] <- values
    }
    mat <- as.matrix(d[species])
    validate(need(!any(is.infinite(mat)),'Infinite abundance values found.'),
             need(!any(mat < 0,na.rm=TRUE),'Negative abundance values found.'))
    list(data=d,species=species,missing=sum(is.na(mat)))
  })
  output$abundance_validation <- renderText({
    if(is.null(input$abundance_file)) return('Upload a site abundance CSV or download the example.')
    x <- abundance_data()
    paste0('Sites: ',nrow(x$data),' | Species columns: ',length(x$species),
           '\\nMissing abundance cells: ',x$missing,
           '\\nZero-abundance sites: ',
           sum(rowSums(x$data[x$species],na.rm=TRUE)==0),
           '\\nData type: site-by-species nonnegative numeric matrix.',
           if(x$missing>0) '\\nChoose an explicit missing-value rule before calculating.' else '')
  })
  output$abundance_preview <- renderTable({
    utils::head(abundance_data()$data,10)
  },striped=TRUE,spacing='xs')
  abundance_computed <- eventReactive(input$run_abundance,{
    x <- abundance_data()
    validate(need(x$missing==0 || identical(input$abundance_missing,'zero'),
      'Missing abundance values detected. Select an explicit rule or correct the CSV.'))
    m <- as.matrix(x$data[x$species])
    if(anyNA(m)) m[is.na(m)] <- 0
    total <- rowSums(m)
    richness <- rowSums(m>0)
    shannon <- numeric(nrow(m))
    simpson <- numeric(nrow(m))
    pielou <- rep(NA_real_,nrow(m))
    for(i in seq_len(nrow(m))) {
      if(total[i]>0) {
        p <- m[i,m[i,]>0]/total[i]
        shannon[i] <- -sum(p*log(p))
        simpson[i] <- 1-sum(p^2)
        if(richness[i]>1) pielou[i] <- shannon[i]/log(richness[i])
      } else {
        shannon[i] <- NA_real_
        simpson[i] <- NA_real_
      }
    }
    data.frame(site_id=x$data$site_id, abundance_total=total,
      richness=richness, shannon=shannon, gini_simpson=simpson,
      effective_species=exp(shannon), pielou_evenness=pielou,
      check.names=FALSE)
  })
  output$abundance_results <- renderTable({
    utils::head(abundance_computed(),100)
  },striped=TRUE,spacing='xs')
  output$abundance_download <- downloadHandler(
    filename=function() 'PlantTraitQuant_abundance_diversity.csv',
    content=function(file) write.csv(abundance_computed(),file,row.names=FALSE,na=''))
  output$pca_example <- downloadHandler(
    filename=function() 'species_traits_example.csv',
    content=function(file) file.copy(example_path,file))
  pca_data <- reactive({
    req(input$pca_file)
    d <- tryCatch(read.csv(input$pca_file$datapath,check.names=FALSE,
                           stringsAsFactors=FALSE,na.strings=c('','NA','NaN'),
                           fileEncoding='UTF-8-BOM'),
                  error=function(e) validate(need(FALSE,paste('Cannot read CSV:',conditionMessage(e)))))
    names(d) <- trimws(names(d))
    validate(need(nrow(d)>0,'CSV has no observations.'),
             need(!anyDuplicated(names(d)),'Duplicate column names detected.'),
             need('species_id' %in% names(d),'CSV must contain species_id.'))
    d$species_id <- trimws(as.character(d$species_id))
    validate(need(!anyNA(d$species_id) && all(nzchar(d$species_id)),'Missing species_id values.'),
             need(!anyDuplicated(d$species_id),'species_id values must be unique.'))
    cols <- setdiff(names(d),'species_id')
    for(nm in cols) {
      if(is.character(d[[nm]])) {
        v <- suppressWarnings(as.numeric(d[[nm]]))
        if(all(is.na(d[[nm]]) | is.finite(v))) d[[nm]] <- v
      }
    }
    numeric_cols <- cols[vapply(d[cols],is.numeric,logical(1))]
    validate(need(length(numeric_cols)>0,'No numeric traits found.'))
    list(data=d,traits=numeric_cols)
  })
  output$pca_trait_choices <- renderUI({
    x <- pca_data()
    selectInput('pca_traits','Select numeric traits (two or more)',
      choices=x$traits,selected=x$traits,multiple=TRUE)
  })
  output$pca_validation <- renderText({
    if(is.null(input$pca_file)) return('Upload a trait CSV or download the example.')
    x <- pca_data()
    chosen <- intersect(input$pca_traits,x$traits)
    if(length(chosen)==0) return(paste('Species:',nrow(x$data),' | Numeric traits:',length(x$traits),
                                     '\\nChoose at least two numeric traits.'))
    good <- complete.cases(x$data[chosen])
    constants <- chosen[vapply(x$data[chosen],function(v) {
      v <- v[is.finite(v)]
      length(v)<2 || length(unique(v))<2
    },logical(1))]
    paste0('Species: ',nrow(x$data),'\\nSelected traits: ',length(chosen),
           '\\nComplete observations: ',sum(good),
           '\\nExcluded rows with missing selected traits: ',sum(!good),
           '\\nConstant/insufficient traits: ',if(length(constants))paste(constants,collapse=', ') else 'none',
           '\\nPCA uses standardized numeric traits and complete cases.')
  })
  pca_computed <- eventReactive(input$run_pca,{
    x <- pca_data()
    chosen <- intersect(input$pca_traits,x$traits)
    validate(need(length(chosen)>=2,'Select at least two numeric traits.'))
    d <- x$data[complete.cases(x$data[chosen]),c('species_id',chosen),drop=FALSE]
    validate(need(nrow(d)>=3,'PCA needs at least three species with complete selected traits.'))
    keep <- chosen[vapply(d[chosen],function(v) all(is.finite(v)) && sd(v)>0,logical(1))]
    validate(need(length(keep)>=2,'At least two nonconstant numeric traits are needed.'))
    fit <- tryCatch(stats::prcomp(d[keep],center=TRUE,scale.=TRUE),
      error=function(e) validate(need(FALSE,paste('PCA failed:',conditionMessage(e)))))
    validate(need(ncol(fit$x)>=2,'At least two PCA axes are required.'))
    scores <- data.frame(species_id=d$species_id,PC1=fit$x[,1],PC2=fit$x[,2],
                         check.names=FALSE)
    loadings <- data.frame(trait=rownames(fit$rotation),
                          PC1=fit$rotation[,1],PC2=fit$rotation[,2],row.names=NULL)
    variance <- data.frame(component=paste0('PC',seq_along(fit$sdev)),
                           variance_percent=100*fit$sdev^2/sum(fit$sdev^2))
    list(scores=scores,loadings=loadings,variance=variance,
         removed=setdiff(chosen,keep))
  })
  output$pca_plot <- renderPlot({
    z <- pca_computed()
    sc <- z$scores
    graphics::par(mar=c(4.5,4.5,2.5,1.5))
    graphics::plot(sc$PC1,sc$PC2,pch=19,col='#176d51',cex=1.1,
                   xlab=paste0('PC1 (',round(z$variance$variance_percent[1],1),'%)'),
                   ylab=paste0('PC2 (',round(z$variance$variance_percent[2],1),'%)'),
                   main='Standardized species trait PCA')
    if(nrow(sc)<=60) graphics::text(sc$PC1,sc$PC2,labels=sc$species_id,pos=3,cex=.65)
    graphics::abline(h=0,v=0,lty=3,col='grey70')
  },width=850,height=550,res=96)
  output$pca_variance <- renderTable(pca_computed()$variance,striped=TRUE)
  output$pca_scores <- renderTable(utils::head(pca_computed()$scores,100),striped=TRUE)
  output$pca_loadings <- renderTable(pca_computed()$loadings,striped=TRUE)
  output$pca_scores_download <- downloadHandler(
    filename=function()'PlantTraitQuant_PCA_scores.csv',
    content=function(file)write.csv(pca_computed()$scores,file,row.names=FALSE))
  output$pca_loadings_download <- downloadHandler(
    filename=function()'PlantTraitQuant_PCA_loadings.csv',
    content=function(file)write.csv(pca_computed()$loadings,file,row.names=FALSE))
  output$download <- downloadHandler(filename=function()'PlantTraitQuant_results.csv',
                                      content=function(file)write.csv(computed(),file,row.names=FALSE,na=''))
}
shinyApp(ui,server)
