# Define server logic required to draw a histogram
server <- function(input, output) {
  # variables for level_5 and metadata
  level_5 <- 0
  metadata <- 0
  
  observeEvent(input$tabs,{
    
    # if the user goes to tab 2
    # functions to accept the csv files and return a data frame for level 5 and metadata
    # displays the tables the user inputs
    # will also generate a rank abundance curve 
    if(input$tabs == "tab2"){
      l5_upload <- reactive({
        inFile1 <- input$file1
        if (is.null(inFile1))
          return(NULL)
        read_csv(inFile1$datapath)
      })
      meta_upload <- reactive({
        inFile2 <- input$file2
        if (is.null(inFile2))
          return(NULL)
        read_csv(inFile2$datapath)
      })
      # upload level 5 and metadata files
      # convert to table format for cards
      observeEvent(input$run, {
        level_5 <- l5_upload()
        level_5 %>% filter_all(any_vars(. != 0))
        metadata <- meta_upload()
        output$l5_contents <- renderTable({head(level_5)})
        output$meta_contents <- renderTable({metadata})
        
        # set column names in metadata to sample and treatment
        colnames(metadata) <- c("sample", "treatment")
        
        sep_taxa <- separate_taxa(level_5, metadata)
        
        # create a global dataset for treatments and samples, indicated by double arrows
        samples <<- sep_taxa$samples
        treatments <<- sep_taxa$treatments
        meta_global <<- metadata
        
        # create global datasets for each taxonomix level
        Phylum <<- create_datasets(sep_taxa, "phylum")
        Class <<- create_datasets(sep_taxa, "class")
        Order <<- create_datasets(sep_taxa, "order")
        Family <<- create_datasets(sep_taxa, "family")
      })
    }
    # server side functions for core taxa visualization
    else if(input$tabs == 'tab3'){
      if(exists('Phylum')){
        
        # creates core dataset based on taxonomic level selected for tab4
        core_data <- reactive ({
          taxa <- input$taxon_core
          my_list <- get(taxa)
          core <- coreTaxa(my_list$ColumnData)
          return(core)
        })
        
        # creates the caption for the core dataset based on the taxonomic level selected
        core_data_caption <- reactive({
          taxa <- input$taxon_core
          my_list <- get(taxa)
          core <- core_taxa(my_list$ColumnData)
          caption <- make_core_caption(my_list$ColumnData,core)
          return(caption)
        })
        # returns caption and core taxa table to output
        output$coreCaption <- renderUI({HTML(core_data_caption())})
        output$coreTaxa <- renderTable({core_data()},digits=0)
      }
    }
    # server side functions for unique taxa visualization
    else if (input$tabs == 'tab4'){
      if(exists('Phylum')){
        
        # create unique taxa output for tab 4 
        # creates unique taxa datset based on taxonomic level selected
        # Create a list with different objects for each treatment
        unique_data <- reactive({
          taxa <- input$taxon_unique
          my_list <- get(taxa)
          unique_taxa_list <- list()
          unique_taxa_list <- uniqueTaxa(my_list$ColumnData, treatments)
          
          return(unique_taxa_list)
        })
        
        # creates a list of the number of output tables based on the number of treatments 
        # tables are labeled with treatment names
        output$unique_tables <-
          renderUI({
            table_output_list <- lapply(treatments$treatment, function(i){
              table_name <- paste0("table", i)
              htmlOutput(tablename)
            })
            # generates a list of HTML tags Shiny uses to categorize each table output
            # allows Shiny to properly render the unique tables
            tagList(table_output_list)
          })
        # render the table
        # use the "observe" command to allow reactive functions to run
        observe({
          unique_data_tables <- unique_data()
          
          if (!length(unique_data_tables))
          {
            output$no_unique_caption <- renderText({"There are no unique taxa."})
          }
          
          for (i in 1:nrow(treatments)) {
            local({
              my_i <- treatments[i,1]
              # create the caption
              num_unique <- length(unique_data_tables[[my_i]]$treatment)
              taxa <- input$taxon_unique
              my_list <- get(taxa)
              total <- nrow(my_list$ColumnData)
              caption_title <- paste('Unique taxa in', my_i, 'treatment:', num_unique, "of", total, "taxa are unique to this treatment")
              
              tablename <- paste0('table', my_i)
              output[[tablename]] <- renderTable(
                {
                  unique_data_tables[[my_i]]
                },
                caption = caption_title, caption.placement = getOption(
                  'xtable.caption.placement', 'top'
                )
              )
            })
          }
        })
      }
      
    }
    # server side functions for rarefaction visualization
    else if (input$tabs == 'tab5'){
      if(exists('Phylum')){
        
        # select data and create the initial rarefaction graph for tab6
        which_initial_rarefaction <- reactive({
          taxa <- input$taxon_rarefy
          bytype <- input$by_type
          my_list <- get(taxa)
          my_data <- my_list$OTU_t
          graphRare(rarecurve(my_data,label=FALSE,step=20), taxa, bytype)
        })
        output$initialRarefaction <- renderPlot({which_initial_rarefaction})
        
        # select data and create even rarefaction graph
        which_even_rarefaction <- reactive ({
          taxa <- input$taxon_rarefy
          bytype <- input$by_type
          my_list <- get(taxa)
          my_data <- my_list$OTU_rarefy_t
          graphRare(rarecurve(my_data,label=FALSE,step=20), taxa, bytype)
        })
        output$evenRarefaction <- renderPlot({which_even_rarefaction})
      }
      
    }
    # server side functions for bar graphs
    else if (input$tabs == 'tab6'){
      if(exists('Phylum')) {
        which_bar_graph <- reactive({
          taxa <- input$taxon_bar
          my_list <- get(taxa)
          if(input$rawrare_bar=='Raw Data'){
            my_data <- my_list$longDataOther
          }
          else {
            my_data <- my_list$longDataRareOther
          }
          if(input$abs_rel == 'Absolute Abundance'){ # absolute abundance
            ggplot(data=my_data, aes(x=Sample, y=Abundance)) +
              geom_bar(aes(fill=my_data[,2]), position='stack', stat='identity')+
              labs(fill=taxa, y='Absolute Abundance')+
              facet_grid(.~treatment, space='free_x', scales='free_x')+
              theme(legend.position='bottom')+
              guides(fill=guide_legend(ncol=2))
          }
          else{
            ggplot(data=my_data, aes(x=Sample, y=Abundance))+
              geom_bar(aes(fill=my_data[,2]), position='fill', stat='identity')+
              labs(fill=taxa, y='Relative Abundance')+
              facet_grid(.~treatment, space='free_x', scales='free_x')+
              theme(legend.position='bottom')+
              guides(fill=guide_legend(ncol=2))
          }
        })
        # determine height of the taxonomy bar graphs in tab3
        bar_graph_height <- reactive({
          taxa <- input$taxa_bar
          my_list <- get(taxa)
          my_data <- my_list$longDataOther
          num_taxa <- length(unique(my_data[,2]))
          height = 300 + num_taxa*10
          return(height)
        })
        # output the taxonomy bar graph 
        observe({output$bargraph <- renderPlot({which_bar_graph()}, height = bar_graph_height())})
      }
    }
    # server side functions for taxa heatmaps
    else if (input$tabs == 'tab7'){
      heat_map_height <- reactive({
        taxa <- input$taxon_heatmap
        my_list <- get(taxa)
        my_data <- my_list$longDataOther
        num_taxa <- nrow(unique(my_data[,1]))
        height <- max(c(225, num_taxa*15))
        return(height)
      })
      heat_map <- reactive({
        taxa <- input$taxon_heatmap
        my_list <- get(taxa)
        if(input$rawrare_heatmap=='Raw Data'){
          my_data <- my_list$longData
          ggplot(my_data, aest(x=Sample, y=.data[[taxa]], fill=Abundance))+
            geom_tile(color='gray')+
            theme(legend.justification='top', axis.text.x=element_text(angle=-90, hjust=0.5))+
            scale_x_discrete(position='top')
        }
      })
      # obeserve the heatmap
      observe({output$heatmap <- renderPlot({heatmap()}, height=heat_map_height())})
    }
    # server side functions for alpha diversity visualization
    else if (input$tabs == 'tab8'){
      # create a boxplot of alpha diversity depending on taxa, data type, and index
      which_alpha_plot <- reactive({
        taxa <- input$taxon_alpha_test
        my_list <- get(taxa)
        if(input$rawrare_alpha == 'Raw Data'){
          my_data <- my_list$DiversityResults
        }
        else{
          my_data <- my_list$DiversityResults_rarefy
        }
        
        which_div <- input$div_measure
        
        data_median <- summarise(group_by(my_data, treatment), MD = round(median(.data[[whichdiv]]),2))
        
        y_label <- names(diversity_choices[grep(which_div, diversity_choices)])
        
        gglot(my_data, aes(x=treatment, y=.data[[which_div]]))+
          geom_boxplot(x='Treatment', y=y_label)+
          geom_text(data=data_median, aes(treatment, MD, label=MD), position=position_dodge(width=0.8), size = 3, vjust=-0.5)
      })
      output$alpha_plots <- renderPlot({which_alpha_plot})
      
      # close alpha diversity observe
    }
    # server side functions for beta diversity visualization
    else if (input$tabs == 'tab9'){
      # return the physeq object
      which_phy_seq <- reactive({
        taxa <- input$taxa_beta_test
        my_list <- get(taxa)
        if(input$rawrare_beta=='Raw Data'){
          my_physeq <- my_list$physeq
        }
        else{
          my_physeq <- my_list$physeq2
        }
        return(my_physeq)
      })
      
      # make ordination data 
      which_ordination_data <- reactive({
        my_ord_data <- ordinate(which_phy_seq(), method=input$ord_method, distance=input$dist_measure)
        return(my_ord_data)
      })
      
      # plot the ordination
      which_ordination_plot <- reactive ({
        if(input$samptreat=='Sample'){
          plot_ordination(which_phy_seq(), which_ordination_data, color = 'treatment') +
            stat_ellipse(type='t')+
            theme_bw()+
            coord_fixed()
        }
      })
      
      # make the ordination caption (NMDS stress)
      which_ordination_caption <- reactive({
        ord_data <- which_ordination_data()
        if(input$ord_method=='NMDS'){
          ord_cap <- paste("NMDS results: stress =", round(ord_data$stress, digits=4))
        }
        # may need to include an else statement here for other ordination methods
      })
      
      # plot ordination in UI
      output$ordination_plot <- renderPlot({which_ordination_plot()})
      output$ordination_caption <- renderText({which_ordination_caption})
      
      # return adonis results
      which_permanova <- reactive({
        taxa <- input$taxon_beta_test
        my_list <- get(taxa)
        if(input$rawrare_beta == 'Raw Data'){
          perm_method <- my_list$OTU_t
        }
        else{
          perm_method <- my_list$OTU_rarefy_t
        }
        
        perm_out <- adonis(perm_method~treatment, data=meta_global, method=input$dist_measure)
        return(as.data.frame(perm_out$aov.tab))
      })
      output$permanova_caption <- renderText({'PERMANOVA results'})
      output$beta_permanova <- renderTable({which_permanova()}, rownames=TRUE)
    }
    # close the observe event
  })
}
