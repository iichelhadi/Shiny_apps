suppressPackageStartupMessages({
  library(tidyverse)
  library(shiny)
  library(shinythemes)
  library(ggfortify)
  library(DESeq2)
  library(pheatmap)
  library(janitor)
  library(umap)
  library(EnhancedVolcano)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(org.Rn.eg.db)
  library(BiocManager)
  library(gtools)
  library(DT)
  library(ggh4x)
  library(ggpubr)
})
# options(repos = BiocManager::repositories())
set.seed(123)

################################################################################
################################## functions ###################################

dds_wide_fn <- function(counts, coldata, group){
  coldata <- coldata %>% mutate_if(is.character,as.factor)
  dds <- DESeqDataSetFromMatrix(counts,coldata,design = ~ group)
  norm <- vst(dds)
  norm_counts <- as.data.frame(assay(norm))
  dim(norm_counts)
  data_wide2 <- norm_counts
  dim(data_wide2)
  data_wide2 <-
    data_wide2 %>% as.data.frame() %>% t %>% cbind(., coldata)
  return(as.data.frame(data_wide2))
}

dds_wide_fn2 <- function(counts, coldata, group){
  coldata <- coldata %>% mutate_if(is.character,as.factor)
  dds <- DESeqDataSetFromMatrix(counts,coldata,design = ~ group)
  norm <- vst(dds)
  norm_counts <- as.data.frame(assay(norm))
  return(norm_counts)
}

dds_long_fn <- function(counts, coldata, group) {
  coldata <- coldata %>% mutate_if(is.character,as.factor)
  dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ group)
  norm <- vst(dds)
  norm_counts <- as.data.frame(assay(norm))
  dim(norm_counts)
  norm_counts$genes <- rownames(norm_counts)
  data2 <- gather(norm_counts, key = 'sample', value = 'expr',-genes)
  data2 <- left_join(data2, coldata, by = 'sample')
  return(as.data.frame(data2))
}

pca_fn <- function(x, group){
  x %>% dplyr::select(where(is.numeric)) %>% 
    dplyr::select(., -caret::nearZeroVar(.)) %>%
    prcomp(scale. = T) %>%
    autoplot(.,data = x, colour = group ,size = 4) +
    theme_bw(base_size = 20)+ guides(color=guide_legend("Cell type"))
}

umap_fn <- function(x, groupvar) {
  num_data_points <- nrow(x)
  n_neighbors <- min(15, num_data_points)
  df <- x %>% dplyr::select(where(is.numeric)) %>% 
    dplyr::select(., -caret::nearZeroVar(.)) %>% scale() %>% 
    umap(.,preserve.seed = T,n_neighbors = n_neighbors)
  df2 <- as.data.frame(df$layout)
  colnames(df2)[c(1,2)] <- c('UMAP_1','UMAP_2')
  df2 <- cbind(df2, x %>% dplyr::select(!where(is.numeric)))
  p <- df2  %>% ggplot(.,aes(x = UMAP_1,y = UMAP_2,color = .data[[groupvar]])) +
    geom_point(size = 4) + labs(x = "UMAP1", y = "UMAP2") +
    theme_bw(base_size = 20) + guides(color = guide_legend(groupvar))
  return(p)
}

expr_fn <- function(x, groupvar, colorvar, gene) {
  x %>% subset(genes %in% c(gene)) %>% 
    ggplot(.,aes(x = .data[[groupvar]], y = expr, fill = .data[[colorvar]])) +
    geom_boxplot(width = 0.4) + geom_jitter(color = 'black',width = .2) +
    scale_x_discrete(guide = "axis_nested") +
    # scale_x_discrete(limits = mixedsort(unique(x$group))) +
    facet_wrap(~ genes ,scales = 'free_y') +
    theme_bw(base_size = 20) + 
    theme(ggh4x.axis.nesttext.x = element_text(angle = 90),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 60,vjust = 1,hjust = 1,size = 10,
                                     face = 'bold'),
          # legend.position = "none",
          strip.text = element_text(size = 15)) +
    xlab(NULL) + ylab('normalized expression')+
    stat_compare_means()
}

diff_res_fn <-  function(counts, coldata, group) {
  dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ group)
  keep <- rowSums(counts(dds,normalized = F) >= 10) >= 3
  dds <- dds[keep,]
  dds <- DESeq(dds)
  return(dds)
}

volcano_fn <- function(res, pval, FC) {
  p <- EnhancedVolcano(
    res,
    lab = rownames(res),
    labSize = 4,
    labFace = "bold",
    x = "log2FoldChange",
    y = "pvalue",
    xlab = bquote( ~ Log[2] ~ "fold change"),
    ylab = bquote( ~ -Log[10] ~ italic(P)),
    pCutoff = pval,
    FCcutoff = FC,
    # xlim = c(-7,5),
    # ylim = c(0,17.5),
    colAlpha  = 0.6,
    legendLabels = c("NS", "Log2 FC", " p-value",
                     " p-value & Log2 FC"),
    legendPosition = "top",
    legendLabSize = 10,
    legendIconSize = 3.0,
    pointSize = 2,
    title = "Differential Expression Genes"
  )
  return(p)
}

heatmap_fn <- function(vsd_sub, tmp_coldata, grouping) {
  p <- pheatmap(
    vsd_sub,
    annotation_col = tmp_coldata[grouping],
    color = colorRampPalette(colors = c('blue', 'white', 'red'))(50),
    cluster_cols = T, show_rownames = T, scale = 'row'
  )
  print(p)
}
################################################################################

################################################################################
################################## UI ##########################################
# Define UI for application
ui <- fluidPage(
  theme = shinytheme('cerulean'),
  navbarPage("RNA-seq analysis app",
    ################################################################################
    ################### tabpanel input Species
             tabPanel(
               'Choose species',
               sidebarPanel(
                 radioButtons(
                   inputId = 'Species',
                   label = 'Species',
                   choices = c(
                     Human = 'org.Hs.eg.db',
                     Mouse = 'org.Mm.eg.db',
                     Rat = 'org.Rn.eg.db'
                   ),
                   selected = 'org.Hs.eg.db'
                 ),
                 tags$hr()
               )),
    ################################################################################
    ################### tabpanel input count data
    tabPanel(
      'Input count data',
      sidebarPanel(
        fileInput(
          inputId = "file1",
          label = "Choose File",
          accept = c("text/csv",
                     "text/comma-separated-values,text/plain",
                     ".csv")
        ),
        radioButtons(
          inputId = 'sep',
          label = 'Separator',
          choices = c(
            Comma = ',',
            Semicolon = ';',
            Tab = '\t',
            Space = ''
          ),
          selected = '\t'
        ),
        tags$hr(),
        checkboxInput("header", "Header", TRUE)
      ),
      mainPanel(h1('Table'),
                tableOutput("contents")),
    ),
    ################################################################################
    ################### tabpanel input coldata
    tabPanel(
      'Input sample data',
      sidebarPanel(
        fileInput(
          inputId = "file2",
          label = "Choose File",
          accept = c("text/csv",
                     "text/comma-separated-values,text/plain",
                     ".csv")
        ),
        radioButtons(
          inputId = 'sep2',
          label = 'Separator',
          choices = c(
            Comma = ',',
            Semicolon = ';',
            Tab = '\t',
            Space = ''
          ),
          selected = '\t'
        ),
        tags$hr(),
        checkboxInput("header2", "Header", TRUE),
        tags$hr(),
        tags$h4('Choose grouping variables'),
        selectInput(
          inputId = 'groupvar1',
          label = 'Grouping variable',
          multiple = F,
          choices = NULL
        )
      ),
      mainPanel(h1('Table'),
                tableOutput("contents2")),
    ),
    ################################################################################
    ################### tabpanel PCA
    tabPanel('PCA plot',
                       sidebarPanel(
                       tags$h4('Choose grouping variables'),
                       selectInput(
                         inputId = 'pcagroup',
                         label = 'Grouping variable',
                         multiple = F,
                         choices = NULL
                       )),
             mainPanel(h1('PCA plot'),
                       plotOutput("pcaplot")),),
    ################################################################################
    ################### tabpanel UMAP
    tabPanel('UMAP plot',
             sidebarPanel(
                       tags$h4('Choose grouping variables'),
                       selectInput(
                         inputId = 'umapgroup',
                         label = 'Grouping variable',
                         multiple = F,
                         choices = NULL
                       )),
             mainPanel(h1('UMAP plot'),
                       plotOutput("umapplot")),),
    ################################################################################
    ######################## volcano plot ##########################################
    ################### tabpanel volcano
    tabPanel(
      'Volcano plot',
      sidebarPanel(
        tags$h4('Choose Fold change cutoff'),
        textInput(
          inputId = "FC",
          label = 'Fold change',
          value = 1
        ),
        tags$h4('Choose p-value cutoff'),
        textInput(
          inputId = "pval",
          label = 'p-value',
          value = 0.05
        )
      ),
      mainPanel(
        h1('Volcano plot'),
        plotOutput("volcanoplot",
                   height = '800px',
                   width = '800px')
      ),
    ),
    
    ################################################################################
    ############################# expr plot #######################################
    ################### tabpanel expr plot
    tabPanel(
      "Expression plot",
      sidebarPanel(
        tags$h4('Choose coloring variable'),
        selectInput(
          inputId = 'colorvar',
          label = 'Coloring variable' ,
          choices = NULL
        ),
        tags$h4('Choose grouping variable'),
        selectInput(
          inputId = 'groupvar',
          label = 'Grouping variable' ,
          choices = NULL
        ),
        tags$h4('Input gene name:'),
        selectizeInput(
          inputId = "gene",
          label = "Select Gene(s)",
          choices = NULL,
          # Populate choices dynamically
          # selected = 'PAX4',
          multiple = TRUE,
          options = list(server = TRUE,
                         search = TRUE)
        ),
      ),
      ## sidebarPanel
      mainPanel(h1('Boxplot'),
                plotOutput("exprplot",
                           height = "800px"),)
    ),
    ################################################################################
    ################### tabpanel Diff expr
    tabPanel(
      "Differential expression",
      sidebarPanel(
        tags$h4('Select test group:'),
        selectInput(
          inputId = "test_group",
          label = NULL,
          selected = NULL,
          selectize = T,
          choices = NULL,
          multiple = F
        ),
        tags$hr(),
        tags$h4('Select reference group:'),
        selectInput(
          inputId = "ref_group",
          label = NULL,
          selected = NULL,
          selectize = T,
          choices = NULL,
          multiple = F
        ),
        actionButton("submitbutton",
                     'Submit', class = "btm btm-primary"),
        tags$hr(),
        downloadButton("dl", "Download"),
      ),
      ## sidebarPanel
      mainPanel(
        h1('Differentially expressed genes'),
        DT::dataTableOutput("contents3"),
      )
    ),
    ################################################################################
    ################### tabpanel heatmap
    tabPanel(
      'Heatmap plot',
      sidebarPanel(
        tags$h4('Choose Fold change cutoff'),
        textInput(
          inputId = "FC2",
          label = 'Fold change',
          value = 3
        ),
        tags$h4('Choose p-value cutoff'),
        textInput(
          inputId = "pval2",
          label = 'p-value',
          value = 0.05
        ),
        tags$h4('Choose grouping variables'),
        selectInput(
          inputId = 'heatlab',
          label = 'Grouping variable',
          multiple = T,
          choices = NULL
        ),
        tags$h4('heatmap dims'),
        textInput(inputId = "heatmap_width", label = "Width (e.g., '400px'):", value = '1000px'),
        textInput(inputId = "heatmap_height", label = "Height (e.g., '400px'):", value = '1000px'),
      ),
      mainPanel(
        h1('heatmap plot'),
        fluidRow(
        div(style = "width: 1000px; height: 1000px; overflow: auto;",
            plotOutput("heatmapplot", width = '100%', height = '100%')
            ) ## div
      ) ## fluidRow
    ) ##mainpanel
      ),## tabpanel
      ################################################################################
  ) ### navbarPage
  ) #### closing UI
  ################################################################################

################################################################################
################################### server #####################################
server = function(input, output, session) {
  ##############################################################################
  ############################ countsInput #####################################
  
  countsInput <- reactive({
    req(input$Species)
    
    file <- input$file1
    if (!is.null(file)) {
      df <- read.table(file$datapath, sep = input$sep,
          header = input$header, row.names = 1) %>% 
        cbind(genes = rownames(.), .)
      
      selected_species <- input$Species
      if (selected_species == 'org.Hs.eg.db') {
        annotation_db <- org.Hs.eg.db
        df$genes <- mapIds(
          x = annotation_db,
          keys = rownames(df),
          column = "SYMBOL",
          keytype = "ENSEMBL",
          multiVals = "first"
        )
      } else if (selected_species == 'org.Mm.eg.db') {
        annotation_db <- org.Mm.eg.db
        df$genes <- mapIds(
          x = annotation_db,
          keys = rownames(df),
          column = "SYMBOL",
          keytype = "ENSEMBL",
          multiVals = "first"
        )
      } else if (selected_species == 'org.Rn.eg.db') {
        annotation_db <- org.Rn.eg.db
        df$genes <- mapIds(
          x = annotation_db,
          keys = rownames(df),
          column = "SYMBOL",
          keytype = "ENSEMBL",
          multiVals = "first"
        )
      }
      df <- df[complete.cases(df$genes),]
      genes <- df$genes
      rownames(df) <- make.names(genes, unique = T)
      df$genes <- NULL
      return(df[1:2000,])
    }
  })
  
  output$contents <- renderTable({
    head(countsInput())
  }, rownames = T)
  
  ##############################################################################
  ############################# coldataInput ###################################
  coldataInput <- reactive({
    file2 <- input$file2
    if(is.null(file2)){ 
      return(NULL) } 
    else {
      coldata_in <- read.table(file2$datapath,sep=input$sep2, header = input$header2,
                               row.names = 1) %>% cbind(sample = rownames(.),.)
    }
  })
  
  output$contents2 <- renderTable({
    coldataInput()
  }, rownames = T)
  
  observe({
    updateSelectInput(session, "ref_group", choices = sort(unique(coldataInput()$group)),
                      selected = sort(unique(coldataInput()$group))[1])
  })
  
  observe({
    updateSelectInput(session, "test_group", choices = sort(unique(coldataInput()$group)),
                      selected = sort(unique(coldataInput()$group))[2])
  })
  
  observe({
    col_names <- names(coldataInput())
    updateSelectInput(session, "groupvar1", choices = col_names)
    updateSelectInput(session, "colorvar", choices = col_names)
    updateSelectInput(session, "groupvar", choices = col_names)
    updateSelectInput(session, "pcagroup", choices = col_names)
    updateSelectInput(session, "umapgroup", choices = col_names)
  })
  
  ##############################################################################
  ################################## PCA #######################################
  
  plot.dat <- reactiveValues(main=NULL, layer1=NULL)
  
  observeEvent(input$pcagroup, {
    req(input$pcagroup)
    req(input$groupvar1)
    plot.dat$main <- reactive({
      pca_fn(dds_wide_fn(countsInput(), coldataInput(), input$groupvar1), input$pcagroup)
    })
  })

  observe({
    print("render PCA")
    output$pcaplot <- renderPlot({ plot.dat$main() })
  })
  
  ##############################################################################
  ################################### UMAP #####################################
  
  plot.dat2 <- reactiveValues(main=NULL, layer1=NULL)

  observeEvent(input$umapgroup, {
    req(input$umapgroup)
    plot.dat2$main <- reactive({
      umap_fn(dds_wide_fn(countsInput(), coldataInput(), input$groupvar1), input$umapgroup)
    })
  })
  
  observe({
    print("render umap")
    output$umapplot <- renderPlot({ plot.dat2$main() })
  })
  
  ##############################################################################
  ############################## diff table ####################################
  
  dds <- reactive({
    if (is.null(input$file1) || is.null(input$file2)) {
      return(NULL)
    } else {
      diff_res_fn(countsInput(),coldataInput(), input$groupvar1)}
  })
  
  res <- reactive({
    req(input$groupvar1)
    req(input$Species)
    req(input$test_group) 
    req(input$ref_group)
    
    diff_res <- results(dds(),alpha = 0.05, contrast = c(input$groupvar1, input$test_group, input$ref_group))
    
    selected_species <- input$Species
    
    if (selected_species == 'org.Hs.eg.db') {
      annotation_db <- org.Hs.eg.db
      diff_res$entrez <- mapIds(
        x = annotation_db,
        keys = diff_res@rownames,
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
      )
    } else if (selected_species == 'org.Mm.eg.db') {
      annotation_db <- org.Mm.eg.db
      diff_res$entrez <- mapIds(
        x = annotation_db,
        keys = diff_res@rownames,
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
      )
    } else if (selected_species == 'org.Rn.eg.db') {
      annotation_db <- org.Rn.eg.db
      diff_res$entrez <- mapIds(
        x = annotation_db,
        keys = diff_res@rownames,
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
      )
    }
    diff_res <- as.data.frame(diff_res[complete.cases(diff_res),])
    return(diff_res)
  })
  
  output$contents3 <- DT::renderDataTable({
    if (input$submitbutton > 0) {
      res()
    }}, rownames = T,)
  output$dl <- downloadHandler(
    filename = function() {"results.csv"},
    content = function(file) {
      write.csv(res(), file)
    }
  )

  ##############################################################################
  ################################ Heatmap plot ################################
  
  plot.dat5 <- reactiveValues(main = NULL, layer1 = NULL)

  observe({
    col_choice <- colnames(coldataInput())
    updateSelectizeInput(
      session,
      inputId = "heatlab",
      choices = colnames(coldataInput()),
      selected = col_choice[1]
    )
  })

  observeEvent(c(input$heatlab,input$pval2,input$FC2), {

    plot.dat5$col_choice <- input$heatlab  # Update the reactive value
    req(input$heatlab)  # Ensure that input$heatlab is not NULL
    req(input$FC2)
    req(input$pval2)
    req(dds())
    res2 <- results(dds(), alpha = 0.05, contrast = c('group',input$test_group,input$ref_group))
    resSig <- subset(res2, (padj < as.numeric(input$pval2) & abs(log2FoldChange) > as.numeric(input$FC2)))
    
    vsd <- dds_wide_fn2(countsInput(), coldataInput(), input$groupvar1)
    vsd_sub <- vsd[rownames(vsd) %in% rownames(resSig),]
    
    plot.dat5$main <-
      function() {
        heatmap_fn(vsd_sub, coldataInput(), input$heatlab)
      }
  })

  observe({
    print("render Heatmap")
    heatmap_width <- as.numeric(gsub("px", "", input$heatmap_width))
    heatmap_height <- as.numeric(gsub("px", "", input$heatmap_height))
    
    output$heatmapplot <- renderPlot({
      plot.dat5$main()
    }, width = heatmap_width, height = heatmap_height)
  })
  
  ##############################################################################
  ################################ Volcano plot ################################
  
  plot.dat4 <- reactiveValues(main=NULL, layer1=NULL)
  
  plot.dat4$main <- reactive({
    volcano_fn(res(), as.numeric(input$pval), as.numeric(input$FC))
  })
  
  observe({
    print("render Volcano")
    output$volcanoplot <- renderPlot({ plot.dat4$main() })
  })
  
  ##############################################################################
  ############################## expression plot ###############################
  # Create a reactive function to store the uploaded count data
  plot.dat3 <- reactiveValues(main = NULL, layer1 = NULL)
  
  plot.dat3$main <- reactive({
    if (!is.null(input$file1) && !is.null(input$file2)) {
      countsInput()
    }
  })

  # Populate the choices in the selectInput widget based on uploaded data
  observeEvent(plot.dat3$main(), {
    gene_choices <- unique(rownames(countsInput()))
    updateSelectizeInput(
      session,
      inputId = "gene",
      choices = gene_choices,
      selected = plot.dat3$selected_gene
    )
  })
  
  ### grouping variable reactive function
  observe({
    col_names <- names(coldataInput())
    updateSelectInput(session, "colorvar", choices = col_names)
    updateSelectInput(session, "groupvar", choices = col_names)
  })
  
  # Create a reactive function to generate the plot based on selected genes
  plot_selected_genes <- reactive({
    req(input$gene)
    req(input$groupvar)
    req(input$colorvar)
    req(input$groupvar1)
    # 
    # if (is.null(plot.dat3$main())) {
    #   return(NULL)
    # }

    expr_fn(
      dds_long_fn(countsInput(), coldataInput(),input$groupvar1),
      input$groupvar,
      input$colorvar,
      input$gene
    )
  })
  
  # Render the plot
  output$exprplot <- renderPlot({
    plot_selected_genes()
  })
  
  ##############################################################################
}

################################################################################
# Run the application 
shinyApp(ui = ui, server = server)

