library(shiny)
library(shinydashboard)
library(magrittr)
source(here::here("plot_theme.R"))
library(ggplot2)
library(shinyWidgets)
library(shinyBS)
library(dplyr)
library(SingleCellExperiment)
shinyOptions(cache = cachem::cache_mem(max_size = 200e6))


# load the sce object.
# reactivity not really important here, because it is used all the time.
sce <- readRDS(here::here("data","filtered_sce_object.rds"))
feature_mapping <- read.csv(here::here("data","features.csv"))

#############
# FUNCTIONS #
#############
annotation_plot_call <- function(data, dimensions){
  ggplot(data, aes(x=!!sym(dimensions[1]), y=!!sym(dimensions[2]), color=celltype_manual))+
    geom_point(shape=16, size=0.42, stroke=0)+
    small_axis(substr(dimensions[1],1,nchar(dimensions[1])-2), fontsize = 8, arrow_length = 20)+
    ggtitle("Control")+
    theme(plot.title = element_blank(),
          legend.text = element_text(size=13.5),
          legend.key.height = unit(0.5, "cm"),
          legend.key.width = unit(0.4,"cm"),
          legend.title = element_blank())+
    scale_y_reverse()+
    scale_color_manual(values=celltype_colors)+
    guides(color = guide_legend(override.aes = list(size = 3)))
}

annotation_plot_fun <- function(srt, dim_red, split_conditions){
  if (dim_red == "PCA"){
    dim_red = "PC"
  }
  
  dimensions <- c(paste0(dim_red, "_1"), paste0(dim_red, "_2"))
  
  # generate the plot
  plot_data <- colData(sce)[, c("celltype_manual", dimensions, "perturbation")] %>% data.frame() %>% 
    dplyr::sample_frac() 
  
  if (split_conditions == "combined"){
    p <- annotation_plot_call(plot_data, dimensions)
  } else if (split_conditions == "ctrl"){
    p <- annotation_plot_call(plot_data %>% filter(perturbation == "ctrl"), dimensions) +
      ggtitle("Control Condition") +
      theme(plot.title = element_text(hjust=0.5, face="bold", size = 15))
  } else if (split_conditions == "notch"){
    p <- annotation_plot_call(plot_data %>% filter(perturbation == "notch"), dimensions) +
      ggtitle("Notch-knockout Condition") +
      theme(plot.title = element_text(hjust=0.5, face="bold", size = 15))
  }
  return(p)
}


feature_plot_call <- function(data,dimensions,gene){
  ggplot(data, aes(x=!!sym(dimensions[1]), y=!!sym(dimensions[2]), color=!!sym(gene)))+
    geom_point(shape=16, size=0.42, stroke=0)+
    small_axis(substr(dimensions[1],1,nchar(dimensions[1])-2), fontsize = 8, arrow_length = 20)+
    scale_y_reverse()+
    theme(legend.key.width = unit(0.9,"line"),
          legend.text = element_text(size=9),
          legend.key.height = unit(0.7,"line"))+
    scale_color_gradient(low="lightgrey", high = "darkred", 
                         name = "Expression\n(logcounts)\n")
}

feature_plot_fun <- function(srt, gene, dim_red, order_cells, split_conditions){
  if (dim_red == "PCA"){
    dim_red = "PC"
  }
  dimensions <- c(paste0(dim_red, "_1"), paste0(dim_red, "_2"))
  mdata <- colData(sce)[, c("celltype_manual", dimensions, "perturbation")] %>% data.frame()
  exp_data <- logcounts(sce)[gene,] %>% data.frame()
  colnames(exp_data) <- gene
  stopifnot(all(rownames(mdata) == rownames(exp_data)))
  
  plot_data <- cbind(mdata, exp_data)
  
  if(order_cells){
    plot_data <- plot_data %>% 
      dplyr::arrange(., !!sym(gene))
  } else {
    plot_data <- plot_data %>% 
      dplyr::sample_frac()
  }
  
  if (split_conditions == "combined"){
    p <- feature_plot_call(plot_data,dimensions,gene)
  } else if (split_conditions == "ctrl"){
    p <- feature_plot_call(plot_data %>% filter(perturbation == "ctrl"),dimensions,gene) +
      ggtitle("Control Condition") +
      theme(plot.title = element_text(hjust=0.5, face="bold", size = 15))
    
  } else if (split_conditions == "notch"){
    p <- feature_plot_call(plot_data %>% filter(perturbation == "notch"),dimensions,gene) +
      ggtitle("Notch-knockout Condition") +
      theme(plot.title = element_text(hjust=0.5, face="bold", size = 15))
  }
  return(p)
}

vln_plot_fun <- function(srt, gene, celltypes, summary_stats, split_conditions){
  mdata <- colData(sce)[, c("celltype_manual", "perturbation")] %>% data.frame()
  exp_data <- logcounts(sce)[gene,] %>% data.frame()
  colnames(exp_data) <- gene
  stopifnot(all(rownames(mdata) == rownames(exp_data)))

  plot_data <- cbind(mdata, exp_data) %>% 
    filter(celltype_manual %in% celltypes)
  
  if (split_conditions == "combined"){
    p <- ggplot(plot_data, aes(x=celltype_manual, y=!!sym(gene), color=celltype_manual))+
      geom_violin(colour = "gray60", alpha = 0.0, 
                  scale = "width", width = 0.8)+
      ggbeeswarm::geom_quasirandom(width = 0.4, alpha = 0.35, 
                                   size=0.4, groupOnX = TRUE, bandwidth = 1)+
      theme_Publication() %+%
      theme(axis.title.x = element_blank(),
            plot.title = element_blank(),
            legend.position = "none")+
      ylab("Expression\n(logcounts)")+
      scale_color_manual(values=celltype_colors)
    
    if (!is.null(summary_stats) & "Include mean" %in% summary_stats){
      p <- p + 
        stat_summary(fun="mean", geom="point", color="black", size=3, shape=8)
    }
    if (!is.null(summary_stats) & "Include median" %in% summary_stats){
      p <- p + 
        stat_summary(fun="median", geom="point", color="black", size=3, shape=2)
    }
  } else if(split_conditions == "split"){
    p <- ggplot(plot_data, aes(x=celltype_manual, y=!!sym(gene), color=perturbation, fill=perturbation))+
      geom_violin(colour = "gray60", alpha = 0.0, 
                  scale = "width", width = 0.65)+
      ggbeeswarm::geom_quasirandom(width = 0.1625, alpha = 0.35, dodge.width=0.65,
                                   size=0.4, groupOnX = TRUE, bandwidth = 1)+
      theme_Publication_side_legend() %+%
      theme(axis.title.x = element_blank(),
            plot.title = element_blank(),
            legend.text = element_text(size=15),
            legend.title = element_text(size=16))+
      ylab("Expression\n(logcounts)")+
      guides(color = guide_legend(title = "Condition", override.aes = list(size = 2,alpha=1, shape=NULL, color=c("#f5756d","#01bdc3"))),
             fill = "none")
    
    if (!is.null(summary_stats) & "Include mean" %in% summary_stats){
      p <- p + 
        stat_summary(fun="mean", geom="point", aes(group = perturbation), color="black", size=3, shape=8, position = position_dodge2(width = 0.65))
    }
    if (!is.null(summary_stats) & "Include median" %in% summary_stats){
      p <- p + 
        stat_summary(fun="median", geom="point", aes(group = perturbation), color="black", size=3, shape=2, position = position_dodge2(width = 0.65))
    }
  }
  return(p)
}

ui <-dashboardPage(
  dashboardHeader(title = "", titleWidth = 100,
                  # Set height of dashboardHeader
                  tags$li(class = "dropdown",
                          tags$style(".main-header {max-height: 20px}"),
                          tags$style(".main-header .logo {height: 20px;}"),
                          tags$style(".sidebar-toggle {height: 20px; padding-top: 1px !important;}"),
                          tags$style(".navbar {min-height:20px !important}")
                  ) 
  ),
  dashboardSidebar(
    tags$style(".left-side, .main-sidebar {padding-top: 20px}"),
    sidebarMenu(
      menuItem("About", tabName = "about", icon = icon("fas fa-circle-info", class = "fa-lg")),
      menuItem("Dataset", tabName = "dataset", icon = icon("fas fa-magnifying-glass", class = "fa-lg")),
      menuItem("Contact", tabName = "contact", icon = icon("far fa-address-book", class = "fa-lg"))
    ), width = 100
  ),
  
  dashboardBody(
    tags$style(
      HTML("
           .col-sm-4 {
             padding: 5px !important; /* Adjust padding value as needed */
           }
           .col-sm-3 {
             padding: 5px !important; /* Adjust padding value as needed */
           }
           .col-sm-8 {
             padding: 5px !important; /* Adjust padding value as needed */
           }
           .box-title {
           
             font-size: 14px !important;
             font-weight: bold !important; /* Add this line to make title bold */

           }
           .box-header {
             padding: 3px 8px !important; /* Adjust padding values as needed */
           }
           .box.box-solid.box-primary {
             margin-bottom: 0px !important; /* Adjust margin value as needed */
           }
           .irs-handle.single {
           width: 5px;
           height: 5px;
           }
           ")
    ),
    tabItems(
      #############
      # ABOUT TAB # 
      #############
      tabItem(tabName = "about",
              h4("Single cell profiling of notch mutant intestinal tumors reveals 
                 Chronophage as a novel transcription factor involved in stem 
                 cell proliferation and differentiation", align="center",
                 style = "width: 80%;  margin: 0 auto; font-weight: bold; font-style: italic;"),
              
              div(style = "margin-bottom: 10px;"),  # Add a margin to create space
              
              h5(HTML("Nick Hirschmüller*, Siamak Redhai*<sup>#</sup>, Shivohum Bahuguna*, 
              Svenja Leible, Michaela Holzem, Tianyu Wang, Fillip Port,
              Lea Bräckow, Wolfgang Huber<sup>#</sup>, Michael Boutros<sup>#</sup>"), 
                 align="center",
                 style = "width: 60%;  margin: 0 auto; color: darkgrey;"),
              
              div(style = "margin-bottom: 5px;"),  # Add a margin to create space
              
              h6("*Contributed equally | #Co-Corresponding", 
                 align="center",
                 style = "width: 80%;  margin: 0 auto; color: darkgrey;"),
              
              div(style = "margin-bottom: 30px;"),  # Add a margin to create space
              
              p("This Shiny App allows users to interactively explore the dataset accompanying the publication mentioned above.", 
                HTML("<br>"), 
                "The source code for the app can be found on ", a("GitHub", href="https://github.com/nickhir/scNotch_shiny"), ". The app is intended to be used in fullscreen on a 1920x1080 monitor.")
              
      ),
      
      
      
      #########
      # PLOTS #
      #########   
      tabItem(tabName = "dataset",
              fluidRow(style = "margin-top: -10px; margin-bottom: -10px; padding: 0px",
                       box(title = "Celltype annotation", 
                           plotOutput("annotation_plot", height = 290), 
                           fluidRow(column(12, align="right", style="margin-top: -295px; margin-right: -200",
                                           downloadButton("download_annotation_plot", label=NULL))),
                           status="primary", solidHeader = TRUE, width = 4),
                       
                       box(title = "Feature Expression", 
                           plotOutput("feature_plot", height = 290), 
                           fluidRow(column(12, align="right", style="margin-top: -295px; margin-right: -200",
                                           downloadButton("download_feature_plot", label=NULL))),
                           status="primary", solidHeader = TRUE, width = 4),
                       
                       
                       box(title = "Settings", 
                           div(style = "margin-bottom: 0px", 
                               selectizeInput("selected_gene", "Gene name", 
                                              choices=c(rownames(sce), feature_mapping$FBid), 
                                              options = list(search_contains = TRUE,  
                                                             dropdownParent = "body",  
                                                             scoreThreshold = 1  
                                              ))), 
                           
                           div(style = "margin-top: 10px", 
                               awesomeRadio("dim_reduction", HTML("Dimensionality reduction"), choices = c("UMAP", "PCA", "PHATE"), selected = "UMAP", inline = TRUE),
                               column(style = "margin-left: -7px", 
                                      width=3, offset = 0,
                                      div(prettyCheckbox("order_cells", "Order cells", TRUE),
                                          prettyCheckbox("split_condition", "Split by condition", FALSE))),
                               column(width=1, offset = 7, 
                                      icon("info-circle", class = "icon-info", id = "order_cells_info"),
                                      div(style = "margin-bottom: 15px;"),  
                                      icon("info-circle", class = "icon-info", id="split_condition_info")
                               )),
                           actionButton("go", "Update All Plots!",
                                        icon = icon("fa-solid fa-gears", class = "fa-lg"),
                                        style = "font-weight: bold; background-color: #14E821; margin-top: 20px;"),
                           
                           status="info", solidHeader = TRUE, width = 3
                       ),
                       bsTooltip(id = "order_cells_info", title = "Should cells be plotted in order of expression?", placement = "right", trigger="hover"),
                       bsTooltip(id = "split_condition_info", title = "Create seperate plot for Ctrl and Notch-knockout condition?", placement = "right", trigger="hover"),
              ),
              
              div(style = "margin-top: 10px",
                  conditionalPanel(
                    condition = "input.split_condition",
                    fluidRow(
                      box(width = 4, title = "Celltype annotation", plotOutput("annotation_plot_split", height = 290), status="primary", solidHeader = TRUE),
                      box(width = 4, title = "Feature Expression", plotOutput("feature_plot_split", height = 290) ,  status="primary", solidHeader = TRUE)
                    )
                  )),
              
              fluidRow(style = "margin-top: -10px; padding: 0px",
                       box(title = "Expression across celltypes", plotOutput("vln_expression", height=200),status="primary", solidHeader = TRUE, width = 8,
                           fluidRow(column(12, align="right", style="margin-top: -205px; margin-right: -200",
                                           downloadButton("download_vln_expression", label=NULL)))),
                       
                       box(title = "Settings", status="info", solidHeader = TRUE, width = 3,
                           div(style = "column-count: 3; margin-bottom: 15px",
                               checkboxGroupInput("celltypes_vln_plot", "Available celltypes",
                                                  choices = c("ISC", "EB", "EEP","dEC", "daEC", "aEC", 
                                                              "mEC", "Copper", "LFC", "pEC", "EE", "MT", "unk"),
                                                  selected = c("ISC", "EB", "EEP","dEC", "daEC", "aEC", 
                                                               "mEC", "Copper", "LFC", "pEC", "EE"))
                           ),
                           
                           div(style="margin-bottom: 0px;",
                               checkboxGroupInput("summary_stats", "Summary statistics",
                                                  choices = c("Include mean", "Include median"),
                                                  selected = NULL, inline = T)
                           )
                       )
              )
      ),
      ###############
      # CONTACT TAB #
      ###############
      tabItem(tabName = "contact",
              h2("Contact Information"),
              p("If you have any questions or found a bug with the app, feel free to reach out to any of the
                people listed below \U1F680."),
              div(style = "margin-bottom: 20px;"),  # Add a margin to create space
              h4("Joint first authors", style="font-weight: bold;"),
              p("Nick Hirschmüller,", 
                a("nick.hirschmueller@embl.de", href = "mailto:nick.hirschmueller@embl.de"), 
                ",",
                a("European Molecular Biology Laboratory (EMBL)", href="https://www.embl.org/"),
                style = "margin-bottom: 3px;"),
              p("Siamak Redhai,", 
                a("siamak.redhai@dkfz-heidelberg.de", href = "mailto:siamak.redhai@dkfz-heidelberg.de"), 
                ",",
                a("German Cancer Research Center (DKFZ)", href="https://www.dkfz.de/en/index.html"),
                style = "margin-bottom: 3px;"),
              p("Shivohum Bahuguna,", 
                a("s.bahuguna@dkfz-heidelberg.de", href = "mailto:s.bahuguna@dkfz-heidelberg.de"), 
                ",",
                a("German Cancer Research Center (DKFZ)", href="https://www.dkfz.de/en/index.html"),
                style = "margin-bottom: 30px;"),
              h4("Corresponding authors", style="font-weight: bold;"),
              p("Siamak Redhai,", 
                a("siamak.redhai@dkfz-heidelberg.de", href = "mailto:siamak.redhai@dkfz-heidelberg.de"), 
                ",",
                a("German Cancer Research Center (DKFZ)", href="https://www.dkfz.de/en/index.html"),
                style = "margin-bottom: 3px;"),
              p("Michael Boutros,", 
                a("m.boutros@dkfz-heidelberg.de", href = "mailto:m.boutros@dkfz-heidelberg.de"), 
                ",",
                a("German Cancer Research Center (DKFZ)", href="https://www.dkfz.de/en/index.html"),
                style = "margin-bottom: 3px;"),
              p("Wolfgang Huber,", 
                a("wolfgang.huber@embl.org", href = "mailto:wolfgang.huber@embl.org"), 
                ",",
                a("European Molecular Biology Laboratory (EMBL)", href="https://www.embl.org/"),
                style = "margin-bottom: 3px;"),
      ) 
    )
  )
)





##########
# SERVER #
##########
server <- function(input, output, session) {
  # server side selection
  updateSelectizeInput(session, 'selected_gene', choices = c(rownames(sce), feature_mapping$FBid), server = TRUE,
                       selected = "esg")
  
  
  # By default, we load some standard plots. So the user sees something upon startup
  output$annotation_plot <- renderPlot(annotation_plot_fun(sce, dim_red = "UMAP",
                                                           split_conditions = "combined")) 
  
  output$feature_plot <- renderPlot(feature_plot_fun(sce, gene="esg",
                                                     dim_red = "UMAP",
                                                     order_cells = TRUE,
                                                     split_conditions = "combined"))
  
  output$vln_expression <- renderPlot(vln_plot_fun(sce, gene="esg",
                                                   celltypes = c("ISC", "EB", "EEP","dEC", "daEC", "aEC", 
                                                                 "mEC", "Copper", "LFC", "pEC", "EE"),
                                                   summary_stats=NULL,
                                                   split_conditions = "combined"))
  
  # generate plots upon button press
  observeEvent(input$go, {
    dim_red <- input$dim_reduction
    order_cells <- input$order_cells
    celltypes <- input$celltypes_vln_plot
    gene_name <- input$selected_gene
    # if the user selected a flybase gene id, we map it to the symbol
    if (grepl("FBgn\\d+", gene_name)){
      gene_name <- feature_mapping %>% filter(FBid == gene_name) %>% pull(symbol)
    }
    summary_stats <- input$summary_stats
    split_condition <- input$split_condition
    
    
    if (!split_condition){
      
      ###################
      # ANNOTATION PLOT #
      ###################
      output$annotation_plot <- renderPlot(annotation_plot_fun(sce, dim_red = dim_red, split_conditions = "combined")) %>%
        bindCache(dim_red, "combined",cache="session")
      
      
      #################
      # Feature PLOT #
      ################
      output$feature_plot <- renderPlot(feature_plot_fun(sce, gene=gene_name,
                                                         dim_red = dim_red,
                                                         order_cells = order_cells,
                                                         split_conditions = "combined"))%>%
        bindCache(dim_red,gene_name, order_cells, "combined", cache="session")
      
      
      
      ################
      # Violin PLOT #
      ###############      
      output$vln_expression <- renderPlot(vln_plot_fun(sce, gene=gene_name,
                                                       celltypes = celltypes,
                                                       summary_stats=summary_stats,
                                                       split_conditions = "combined"))%>%
        bindCache(celltypes,gene_name, summary_stats, "combined", cache="session")
      
    }
    
    # if split condition, generate different plots
    if (split_condition){
      
      ###################
      # ANNOTATION PLOT #
      ###################
      output$annotation_plot <- renderPlot(annotation_plot_fun(sce, dim_red = dim_red, split_conditions = "ctrl")) %>%
        bindCache(dim_red, "ctrl",cache="session")
      
      output$annotation_plot_split <- renderPlot(annotation_plot_fun(sce, dim_red = dim_red, split_conditions = "notch")) %>%
        bindCache(dim_red, "notch",cache="session")
      
      #################
      # Feature PLOT #
      ################
      output$feature_plot <- renderPlot(feature_plot_fun(sce, gene=gene_name,
                                                         dim_red = dim_red,
                                                         order_cells = order_cells,
                                                         split_conditions = "ctrl")) %>%
        bindCache(dim_red,gene_name, order_cells, "ctrl", cache="session")
      
      output$feature_plot_split <- renderPlot(feature_plot_fun(sce, gene=gene_name,
                                                               dim_red = dim_red,
                                                               order_cells = order_cells,
                                                               split_conditions = "notch")) %>%
        bindCache(dim_red,gene_name, order_cells, "notch", cache="session")
      
      ################
      # Violin PLOT #
      ###############    
      output$vln_expression <- renderPlot(vln_plot_fun(sce, gene=gene_name,
                                                       celltypes = celltypes,
                                                       summary_stats=summary_stats,
                                                       split_conditions = "split")) %>%
        bindCache(celltypes,gene_name, summary_stats, "split", cache="session")
      
    }
  })
  
  
  
  
  
  
  #################
  # DownloadCalls #
  #################
  output$download_annotation_plot <- downloadHandler(
    filename = function() {
      paste0("annotation_", input$dim_reduction, ".png")
    },
    
    content = function(file) {
      png(file, width = 4.3, height = 3.3, res = 550, units = "in")
      print(annotation_plot_fun(sce, dim_red = input$dim_reduction, 
                                split_conditions = ifelse(input$split_condition, "ctrl","combined")))
      dev.off()
    }
  )
  
  output$download_feature_plot <- downloadHandler(
    filename = function() {
      paste0("feature_expression_", input$selected_gene, ".png")
    },
    
    content = function(file) {
      png(file, width = 4.3, height = 3.3, res = 550, units = "in")
      print(feature_plot_fun(sce, gene=input$selected_gene,
                             dim_red = input$dim_reduction,
                             order_cells = input$order_cells, 
                             split_conditions = ifelse(input$split_condition, "ctrl","combined")))
      dev.off()
    }
  )
  
  output$download_vln_expression <- downloadHandler(
    filename = function() {
      paste0("celltype_expression_", input$selected_gene, ".png")
    },
    
    content = function(file) {
      png(file, width = 7, height = 2, res = 550, units = "in")
      print(vln_plot_fun(sce, gene=input$selected_gene,
                         celltypes = input$celltypes_vln_plot,
                         summary_stats=input$summary_stats, 
                         split_conditions = ifelse(input$split_condition, "ctrl","combined")))
      dev.off()
    }
  )
}

shinyApp(ui, server) 