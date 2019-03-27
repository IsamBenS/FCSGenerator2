server <- function(input, output, session)
{
    useShinyjs()
    #======================================================================================================================
    #======================REACTIVE VALUES=================================================================================
    #======================================================================================================================
    
    app.variables <- reactiveValues(
        fcs.files = NULL,
        fcs.file.info.table = NULL,
        sets.list = NULL,
        populations.list = NULL,
        output.matrices = NULL,
        reduction.percentages = NULL,
        loaded.mutant = NULL,
        used.events = NULL,
        log.text = NULL
    )
    
    env.var <- reactiveValues(
        exp.folder = "/media/data/html/INPUT/",
        tool.wd = system.file("shinyApp", "app", package = "FCSGenerator2")
    )
    
    
    #======================================================================================================================
    #============================FUNCTIONS=================================================================================
    #======================================================================================================================
    
    maj.sets.list <- function() #MAJ DES SETS
    {
        if(length(app.variables$fcs.files)>0)
        {
            lapply(1:length(app.variables$fcs.files), function(i)
            {
                if(length(app.variables$fcs.files[[i]]))
                {
                    col.id <- which(unlist(app.variables$fcs.files[[i]][["markers"]]==
                                               input[[paste0("t_1_pop_sel_",i)]]))
                    if(length(col.id)>0)
                    {
                        updateSelectInput(session, paste0("t_1_set_sel_",i), "Change Cohort", 
                                          choices = app.variables$sets.list,
                                          selected = app.variables$fcs.files[[i]][["set"]])
                    }
                }
            })
        }
    }
    
    update.sets.list <- function(tab.id)
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs==paste0("t_",tab.id))
        {
            selected.set <- app.variables$sets.list
            if(!is.na(input[[paste0("t_",tab.id,"_set_list")]]) && input[[paste0("t_",tab.id,"_set_list")]]!="")
            {
                selected.set <- input[[paste0("t_",tab.id,"_set_list")]]
            }
            updateSelectInput(session, paste0("t_",tab.id,"_set_list"), "Select a Cohort", 
                              choices=app.variables$sets.list,
                              selected=selected.set)
        }
    }
    
    update.files.list <- function(tab.id)
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs==paste0("t_",tab.id)) 
        {
            files.list <- list()
            selected.file <- list()
            if(length(input[[paste0("t_",tab.id,"_set_list")]])>0)
            {
                for(i in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[i]])>0)
                    {
                        if(app.variables$fcs.files[[i]][["set"]] ==  input[[paste0("t_",tab.id,"_set_list")]])
                        {
                            files.list[[length(files.list)+1]] <- names(app.variables$fcs.files)[[i]]
                        }
                    }
                }
            }
            selected.file <- files.list
            if(!is.na(input[[paste0("t_",tab.id,"_file_list")]]) && input[[paste0("t_",tab.id,"_file_list")]]!="")
            {
                selected.file <- as.list(unlist(input[[paste0("t_",tab.id,"_file_list")]])[input[[paste0("t_",tab.id,"_file_list")]]%in%files.list])
            }
            updateSelectInput(session, paste0("t_",tab.id,"_file_list"), "Select a File", choices=files.list, selected=selected.file)
        }
    }
    
    init.file <- function(x, tmp.name, mark.list, pop.names = NULL, compensated=F, transformed=F, group = "SRC")
    {
        if( is.null(app.variables$fcs.files) )    
        {
            app.variables$fcs.files <<- list()
            app.variables$fcs.files.backup <<- list()
            app.variables$sets.list <<- list()
            app.variables$populations.list <<- list()
        }
        markers.default.values <- lapply(1:length(mark.list), function(i)
        {
            return(list(mean(x@exprs[,i]), sd(x@exprs[,i])))
        })
        pop.col <- ncol(x@exprs)
        pop.list <- sort(as.numeric(unique(x@exprs[,pop.col])))
        pop.events.ids <- list()
        if(length(pop.list)<100)
        {
            pop.events.ids <- lapply(1:length(pop.list), function(i)
            {
                return(unlist(which(x@exprs[,pop.col]==pop.list[[i]])))
            })
            if(is.null(pop.names))
            {
                names(pop.events.ids) <- 1:length(pop.list)
            }
            else
            {
                names(pop.events.ids) <- pop.names
            }
        }
        
        
        app.variables$fcs.files[[tmp.name]] <<- list()
        app.variables$fcs.files[[tmp.name]][["file"]] <<- x
        app.variables$fcs.files[[tmp.name]][["markers"]] <<- mark.list
        app.variables$fcs.files[[tmp.name]][["name"]] <<- tmp.name
        app.variables$fcs.files[[tmp.name]][["type"]] <<- group
        app.variables$fcs.files[[tmp.name]][["mutant_ui"]] <<- TRUE
        app.variables$fcs.files[[tmp.name]][["set"]] <<- tmp.name
        app.variables$fcs.files[[tmp.name]][["markers_default_values"]] <<- list()
        app.variables$fcs.files[[tmp.name]][["populations_column"]] <<- ncol(x@exprs)
        app.variables$fcs.files[[tmp.name]][["source_ctrl"]] <<- tmp.name
        app.variables$fcs.files[[tmp.name]][["comp"]] <<- compensated
        app.variables$fcs.files[[tmp.name]][["transf"]] <<- transformed
        #==
        app.variables$output.matrices[[tmp.name]] <<- x@exprs
        #==
        app.variables$used.events[[tmp.name]] <<- rep(T,nrow(x@exprs))
        #==
        app.variables$sets.list[[length(app.variables$sets.list)+1]] <<- tmp.name
        #==
        app.variables$populations.list[[tmp.name]] <<- pop.events.ids
        #==
        app.variables$reduction.percentages[[tmp.name]] <<- list()
        
        # file.vec <- matrix(ncol=5,nrow=1)
        # file.vec[1,1] <- paste0(basename(substr(f,1,nchar(f)-4)), "_", length(app.variables$fcs.files)-1)
        # file.vec[1,2] <- trunc(file.size(f)/1024/1024*1000)/1000
        # file.vec[1,3] <- ncol(x)
        # file.vec[1,4] <- nrow(x)
        # file.vec[1,5] <- nrow(x)
        # 
        # app.variables$file.info.table <<- rbind(app.variables$file.info.table,
        #                                         file.vec)
    }
    
    update.log <- function(text)
    {
        app.variables$log.text <- paste(app.variables$log.text, text, sep = "\n")
    }
    
    list.to.string <- function(list)
    {
        new.string <- ""
        for(i in 1:length(list))
        {
            new.string <- paste(new.string, list[[i]], sep=", ")
        }
        return(new.string)
    }
    
    
    #======================================================================================================================
    #===============================CODE===================================================================================
    #======================================================================================================================
    
    #=======================================LOAD FILES===================================================
    
    observeEvent(input$t_1_select,  #SELECT FILES
    {
        shinyjs::disable("t_1_select")
        m <- matrix(nrow=1,ncol=2)
        m[1,1] = "FlowFrames"
        m[1,2] = "*.csv;*.fcs"
        #==
        temp.files <- choose.files(filters = m,multi = T)
        
        if(length(temp.files) > 0)
        {
            progress <- Progress$new()
            progress$set(message = "Loading Files", value = 0)
            update.log("LOADING Files")
            lapply(temp.files, function(f)
            {
                l <- length(f)
                x <- NULL
                mark.list <- list()
                if(grepl("csv",f))
                {
                    x <- as.matrix(read.csv(f))
                    x <- flowFrame(x)
                    lapply(1:ncol(x@exprs), function(i)
                    {
                        nx <- x@description[[paste0("$P",i,"S")]]
                        if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                        {
                            mark.list[[i]] <<- nx
                        }
                        else
                        {
                            mark.list[[i]] <<- colnames(x)[i]
                        }
                    })
                }
                else
                {
                    x <- read.FCS(f,emptyValue = FALSE)
                    lapply(1:ncol(x@exprs), function(i)
                    {
                        nx <- x@description[[paste0("$P",i,"S")]]
                        if(!is.null(nx) && !is.na(nx) && nx != "" && nx != " ")
                        {
                            mark.list[[i]] <<- nx
                        }
                        else
                        {
                            mark.list[[i]] <<- colnames(x)[i]
                        }
                    })
                }
                
                tmp.name <- paste0(basename(substr(f,1,nchar(f)-4)),"_",length(app.variables$fcs.files))
                init.file(x, tmp.name, mark.list)
                progress$inc(1/length(temp.files), detail = paste(f, "loaded"))
                update.log(paste("========",f, "loaded"))
            })
            progress$set(message = "Files loaded", value = 1)
            update.log("======== FILES LOADED")
            update.log("")
            delay(700, progress$close())
        }
        else
        {
            showNotification("NO FILES SELECTED", duration=5, type="error")
        }
        
        maj.sets.list()
        delay(500, shinyjs::enable("t_1_select"))
    })
    
    observeEvent(input$t_1_create,  #CREATE FILES
    {
        shinyjs::disable("t_1_create")
        if(length(input$t_1_nmb_files)>0 && !is.na(input$t_1_nmb_files) && 
           as.numeric(input$t_1_nmb_files)>0)
        {
            if(length(input$t_1_nmb_populations)>0 && !is.na(input$t_1_nmb_populations) && 
               as.numeric(input$t_1_nmb_populations)>0)
            {
                if(length(input$t_1_nmb_markers)>0 && !is.na(input$t_1_nmb_markers) && 
                   as.numeric(input$t_1_nmb_markers)>0)
                {
                    if(length(input$t_1_nmb_events)>0 && !is.na(input$t_1_nmb_events) && 
                       as.numeric(input$t_1_nmb_events)>0)
                    {
                        if(length(input$t_1_nmb_rare_populations)>0 && !is.na(input$t_1_nmb_rare_populations) && 
                           as.numeric(input$t_1_nmb_rare_populations)>=0)
                        {
                            if(length(input$t_1_freq_rare_populations)>0 && !is.na(input$t_1_freq_rare_populations) && 
                               as.numeric(input$t_1_freq_rare_populations)>=0)
                            {
                                progress <- Progress$new()
                                progress$set(message = "Generating Files", value = 0)
                                update.log("GENERATING FILES")
                                
                                nmb.files <- input$t_1_nmb_files
                                nmb.events <- input$t_1_nmb_events
                                nmb.markers <- input$t_1_nmb_markers
                                nmb.populations <- input$t_1_nmb_populations
                                nmb.rare.populations <- input$t_1_nmb_rare_populations
                                rare.pop.freq <- input$t_1_freq_rare_populations
                                pop.names <- lapply(1:(nmb.populations-nmb.rare.populations), function(pop)
                                {
                                    return(input[[paste0("t_1_pop_",pop,"_name")]])
                                })
                                tmp.freq <- c()
                                if(nmb.rare.populations>0)
                                {
                                    pop.names <- list(unlist(pop.names), 
                                                      paste0("MIN_POP_", 1:nmb.rare.populations))
                                    
                                    tmp.freq <- c(0, sort(sample(2:(as.integer(nmb.events*rare.pop.freq/100)-1), nmb.rare.populations-1)),
                                                  as.integer(nmb.events*rare.pop.freq/100))
                                }
                                
                                min.freq.list <- sapply(1:nmb.populations, function(i)
                                {
                                    if(i<=(nmb.populations-nmb.rare.populations))
                                    {
                                        return(as.numeric(input[[paste0("t_1_pop_",i,"min__freq")]]))
                                    }
                                    else
                                    {
                                        tmp.val <- (tmp.freq[i-nmb.populations+nmb.rare.populations+1] - 
                                            tmp.freq[i-nmb.populations+nmb.rare.populations]+1)/nmb.events
                                        return(tmp.val)
                                    }
                                })
                                
                                max.freq.list <- sapply(1:nmb.populations, function(i)
                                {
                                    if(i<=(nmb.populations-nmb.rare.populations))
                                    {
                                        return(as.numeric(input[[paste0("t_1_pop_",i,"max__freq")]]))
                                    }
                                    else
                                    {
                                        tmp.val <- (tmp.freq[i-nmb.populations+nmb.rare.populations+1] - 
                                                        tmp.freq[i-nmb.populations+nmb.rare.populations]+1)/nmb.events*100
                                        return(tmp.val)
                                    }
                                })
                                
                                for(current.file in 1:nmb.files)
                                {
                                    tmp.name <- paste0("GEN_",current.file,"__",nmb.events,"_",
                                                       nmb.markers,"_",nmb.populations,"_",nmb.rare.populations,"_",
                                                       trunc(as.numeric(Sys.time())))
                                    x <- generate.fcs(nmb.events,nmb.markers,nmb.populations,nmb.rare.populations,rare.pop.freq,
                                                      min.freq.list,max.freq.list)
                                    if(!is.null(x))
                                    {
                                        mark.list <- colnames(x@exprs)
                                        pop.names <- unlist(pop.names)
                                        init.file(x, tmp.name, mark.list, pop.names = pop.names, compensated = T, transformed = T)
                                    }
                                    progress$inc(1/nmb.files, detail = paste(tmp.name, "generated"))
                                    update.log(paste("========", tmp.name, "generated"))
                                }
                                progress$set(message = "Files generated", value = 1)
                                update.log("======== FILES GENERATED")
                                update.log("")
                                delay(700, progress$close())
                            }
                        }
                    }
                }
            }
        }
        delay(500, shinyjs::enable("t_1_create"))
    })
    
    # file.table.fct <- function() #FILES INFORMATION TABLE - UPDATE
    # {
    #     tmp.mat <- app.variables$file.info.table
    #     l <- ncol(tmp.mat)
    #     if(length(l)>0 && l==1)
    #     {
    #         tmp.mat <- t(tmp.mat)
    #         colnames(app.variables$file.info.table) <<- c("Filename", "Size (Mo)", "Number of markers", "Number of events", "Subsample", "TYPE")
    #     }
    #     
    #     return(tmp.mat)
    # }
    # 
    # output$t_1_fileInfo <- renderTable(file.table.fct()) #FILES INFORMATION
    
    output$t_1_files_main <- renderUI(  #CREATION DE L'UI DE CHAQUE FICHIER ET OBSERVE EVENT SUR LA CREATION DES SETS
    {
        file.ui <- NULL
        if(length(app.variables$fcs.files)>0 && input$tabs=="t_1")
        {
            file.ui <- lapply(1:length(app.variables$fcs.files), function(i)
            {
                tmp.ui <- NULL
                if(length(app.variables$fcs.files[[i]]) > 0 )
                {
                    file.type <- app.variables$fcs.files[[i]][["type"]]
                    file.type <- ifelse(file.type=="MUT", 2, 1)
                    #==
                    pop.col <- app.variables$fcs.files[[i]][["populations_column"]]
                    pop.col <- app.variables$fcs.files[[i]][["markers"]][[pop.col]]
                    if(length(input[[paste0("t_1_pop_sel_",i)]])>0 && input[[paste0("t_1_pop_sel_",i)]]!="" && 
                       input[[paste0("t_1_pop_sel_",i)]]%in%app.variables$fcs.files[[i]][["markers"]])
                    {
                        pop.col <- input[[paste0("t_1_pop_sel_",i)]]
                    }
                    #==
                    visualized.markers <- app.variables$fcs.files[[i]][["markers"]]
                    if(!is.na(input[[paste0("t_1_viewed_markers_",i)]]) && length(input[[paste0("t_1_viewed_markers_",i)]])>0)
                    {
                        visualized.markers <- 
                            as.list(unlist(input[[paste0("t_1_viewed_markers_",i)]])[which(input[[paste0("t_1_viewed_markers_",i)]]%in%
                                                                                               app.variables$fcs.files[[i]][["markers"]])])
                    }
                    visualized.markers <- as.list(unlist(visualized.markers)[which(visualized.markers!=pop.col)])
                    #==
                    tmp.ui <- tagList(
                        box(
                            width = 8, collapsible=T, style="min-height:15vh",
                            title = app.variables$fcs.files[[i]][["name"]],
                            style="padding-left:5%",
                            plotlyOutput(paste0("t_1_file_",i), width = "100%")
                        ),
                        box(
                            width = 2, style="height:auto",
                            selectInput(paste0("t_1_pop_sel_",i), "Populations Column",
                                        choices=app.variables$fcs.files[[i]][["markers"]],
                                        selected = pop.col),
                            selectInput(paste0("t_1_cm_sel_",i), "Control/Mutant", 
                                        choices=list("Control"=1,"Mutant"=2),
                                        selected = file.type),
                            selectInput(paste0("t_1_viewed_markers_",i), "Visualized Markers", 
                                        choices=app.variables$fcs.files[[i]][["markers"]], 
                                        selected=visualized.markers, multiple = T),
                            actionButton(paste0("t_1_file_remove_",i), "Remove File", width="60%", style="margin-left:20%")
                        ),
                        box(
                            width = 2, style="auto",
                            # textInput(paste0("t_1_set_cr_",i), "New Cohort", width="90%"),
                            # actionButton(paste0("t_1_set_add_",i), "Create", width="90%"),
                            selectInput(paste0("t_1_set_sel_",i), "Change Cohort", 
                                        choices=app.variables$sets.list,
                                        selected = app.variables$fcs.files[[i]][["set"]])
                        )
                    )
                    
                    # observeEvent(input[[paste0("t_1_set_add_",i)]],
                    # {
                    #      if(length(input[[paste0("t_1_set_cr_",i)]])>0)
                    #      {
                    #          tmp.set <- input[[paste0("t_1_set_cr_",i)]]
                    #          
                    #          if(!(tmp.set %in% unlist(app.variables$sets.list)))
                    #          {
                    #              app.variables$sets.list[[length(app.variables$sets.list)+1]] <<- input[[paste0("t_1_set_cr_",i)]]
                    #              app.variables$fcs.files[[i]][["set"]] <<- input[[paste0("t_1_set_cr_",i)]]
                    #              maj.sets.list()
                    #          }
                    #      }
                    #  }, once = T)
                    
                    observeEvent(input[[paste0("t_1_file_remove_",i)]],
                    {
                        update.log(paste(app.variables$fcs.files[[i]][["name"]], "removed"))
                        update.log("")
                        app.variables$fcs.files[[i]] <<- list()
                        app.variables$output.matrices[[i]] <<- list()
                        app.variables$populations.list[[i]] <<- list()
                        app.variables$reduction.percentages[[i]] <<- list()
                        removeUI(paste0("#t_1_file_",i,"_ui"))
                    }, once = T)
                    
                    return(tmp.ui)
                }
            })
        }
        return(file.ui)
    })
    
    observe( #MAJ DES HEATMAPS DES POPS & LISTE MEAN/SD PAR POP
    {
        if(length(app.variables$fcs.files)>0 && input$tabs=="t_1")
        {
            lapply(1:length(app.variables$fcs.files), function(i)
            {
                if(length(app.variables$fcs.files[[i]]) > 0 && !is.na(input[[paste0("t_1_viewed_markers_",i)]]) && 
                   length(input[[paste0("t_1_viewed_markers_",i)]])>0)
                {
                    pop.ids <- app.variables$populations.list[[i]]
                    exp.mat <- app.variables$output.matrices[[i]]
                    col.id <- app.variables$fcs.files[[i]][["populations_column"]]
                    viewed.markers <- as.numeric(unlist(which(app.variables$fcs.files[[i]][["markers"]]%in%
                                                                  input[[paste0("t_1_viewed_markers_",i)]])))
                    viewed.markers <- viewed.markers[viewed.markers!=col.id]
                    
                    file.plot <- NULL
                    tmp.val <- create.pop.table.from.fcs.matrix(exp.mat, pop.ids, viewed.markers, col.id)
                    if(!is.null(tmp.val[[1]]))
                    {
                        file.plot <- heatmaply(tmp.val[[2]], Rowv = T, Colv="Rowv", dendrogram = "none")
                    }
                    output[[paste0("t_1_file_",i)]] <- renderPlotly(file.plot)
                }
            })
        }
    })
    
    observe( #GESTION DES SETS
    {
        if(length(app.variables$fcs.files)>0 && input$tabs=="t_1")
        {
            lapply(1:length(app.variables$fcs.files), function(i)
            {
                if(length(app.variables$fcs.files[[i]])>0)
                {
                    if(length(app.variables$fcs.files[[i]])>0 && length(input[[paste0("t_1_set_sel_",i)]])>0)
                    {
                        if(app.variables$fcs.files[[i]][["set"]]!=input[[paste0("t_1_set_sel_",i)]])
                        {
                            update.log(paste0(app.variables$fcs.files[[i]][["name"]], ": cohort changed"))
                            update.log(paste("======== from", app.variables$fcs.files[[i]][["set"]]))
                            update.log(paste("======== to", input[[paste0("t_1_set_sel_",i)]]))
                            update.log("")
                            #==
                            app.variables$fcs.files[[i]][["set"]] <<- input[[paste0("t_1_set_sel_",i)]]
                            app.variables$fcs.files[[i]][["source_ctrl"]] <<- input[[paste0("t_1_set_sel_",i)]]
                        }
                    }
                }
            })
        }
    })
    
    observe( #UPDATE POP COLUMN
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_1")
        {
            lapply(1:length(app.variables$fcs.files), function(f.id)
            {
                if(length(app.variables$fcs.files[[f.id]])>0)
                {
                    if(length(input[[paste0("t_1_pop_sel_",f.id)]])>0 && input[[paste0("t_1_pop_sel_",f.id)]] != "" && !is.na(input[[paste0("t_1_pop_sel_",f.id)]]))
                    {
                        pop.col <- as.numeric(which(app.variables$fcs.files[[f.id]][["markers"]]==input[[paste0("t_1_pop_sel_",f.id)]])[[1]])
                        #==
                        if(!is.na(pop.col) && as.numeric(pop.col)>0 && pop.col != app.variables$fcs.files[[f.id]][["populations_column"]])
                        {
                            pop.list <- as.numeric(unique(app.variables$output.matrices[[f.id]][,pop.col]))
                            if(length(pop.list)<=100)
                            {
                                app.variables$fcs.files[[f.id]][["populations_column"]] <<- pop.col
                                pop.events.ids <- lapply(1:length(pop.list), function(i)
                                {
                                    return(unlist(which(app.variables$output.matrices[[f.id]][,pop.col]==pop.list[[i]])))
                                })
                                names(pop.events.ids) <- 1:length(pop.list)
                                app.variables$populations.list[[f.id]] <<- pop.events.ids
                            }
                        }
                    }
                }
            })
        }
    })
    
    observe( #UPDATE FILE TYPE
    {
        if(length(app.variables$fcs.files)>0 && input$tabs=="t_1")
        {
            lapply(1:length(app.variables$fcs.files), function(i)
            {
                if(length(app.variables$fcs.files[[i]]) > 0 && length(input[[paste0("t_1_cm_sel_",i)]])>0 && input[[paste0("t_1_cm_sel_",i)]]!="")
                {
                    file.type <- input[[paste0("t_1_cm_sel_",i)]]
                    if(file.type==2)
                    {
                        app.variables$fcs.files[[i]][["type"]] <<- "MUT"
                    }
                    else
                    {
                        app.variables$fcs.files[[i]][["type"]] <<- "CTRL"
                    }
                    
                }
            })
        }
    })
    
    output$t_1_pop_list <- renderUI(
    {
        pop.ui <- list()
        if(length(input$t_1_nmb_populations)>0 && !is.na(input$t_1_nmb_populations) && 
           as.numeric(input$t_1_nmb_populations)>0 && input$t_1_tb=="A" &&
           length(input$t_1_nmb_rare_populations)>0 && !is.na(input$t_1_nmb_rare_populations) && 
           as.numeric(input$t_1_nmb_rare_populations)>=0 && input$tabs=="t_1")
        {
            pop.ui <- lapply(1:(as.numeric(input$t_1_nmb_populations)-as.numeric(input$t_1_nmb_rare_populations)), function(i)
            {
                nmb.non.rare.pop <- as.numeric(input$t_1_nmb_populations)-as.numeric(input$t_1_nmb_rare_populations)
                effective.perc <- 100
                if(!is.na(input$t_1_freq_rare_populations))
                {
                    effective.perc <- max(0,100-input$t_1_freq_rare_populations)
                }
                val <- tagList(
                    column(
                        width=4,
                        textInput(paste0("t_1_pop_",i,"_name"), "Population Name", value = i)
                    ),
                    column(
                        width=3,
                        numericInput(paste0("t_1_pop_",i,"min__freq"), "Min Frequency", value = trunc(effective.perc/nmb.non.rare.pop*100)/100)
                    ),
                    column(
                        width=3,
                        numericInput(paste0("t_1_pop_",i,"max__freq"), "Max Frequency", value = trunc(effective.perc/nmb.non.rare.pop*100)/100+1)
                    )
                )
                return(val)
            })
            return(pop.ui)
        }
        return(pop.ui)
    })
    
    observe( #UPDATE NUMBER MIN MARKERS
    {
        if(!is.na(input$t_1_nmb_populations) && !is.na(input$t_1_nmb_markers) && input$tabs=="t_1")
        {
            if(as.numeric(input$t_1_nmb_populations) > 2^as.numeric(input$t_1_nmb_markers))
            {
                updateNumericInput(session, "t_1_nmb_markers", "Number of Markers", value=as.integer(log2(as.numeric(input$t_1_nmb_populations)))+1)
            }
        }
    })
    
    
    
    
    
    
    
    #=======================================MODIFY POPULATIONS===========================================
    
    
    
    
    
    
    
    observe( #CLEAR UI
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_2")
        {
            if(length(input[["t_2_file_list"]])>0 && input[["t_2_file_list"]] != "")
            {
                f.name <- input[["t_2_file_list"]]
                if(length(app.variables$fcs.files[[f.name]])>0)
                {
                    shinyjs::show("t_2_plot_div")
                    shinyjs::show("t_2_right_col")
                    shinyjs::enable("t_2_pop_inc")
                    shinyjs::enable("t_2_pop_red")
                    shinyjs::enable("t_2_pop_new")
                }
                else
                {
                    hide("t_2_plot_div")
                    hide("t_2_right_col")
                    shinyjs::disable("t_2_pop_inc")
                    shinyjs::disable("t_2_pop_red")
                    shinyjs::disable("t_2_pop_new")
                }
            }
            else
            {
                hide("t_2_plot_div")
                hide("t_2_right_col")
                shinyjs::disable("t_2_pop_inc")
                shinyjs::disable("t_2_pop_red")
                shinyjs::disable("t_2_pop_new")
            }
        }
        else
        {
            hide("t_2_plot_div")
            hide("t_2_right_col")
            shinyjs::disable("t_2_pop_inc")
            shinyjs::disable("t_2_pop_red")
            shinyjs::disable("t_2_pop_new")
        }
    })
    
    observe( #UPDATE SETS LIST
    {
        update.sets.list(2)
    })
    
    observe( #UPDATE FILES LIST
    {
        update.files.list(2)
    })
    
    observe( #UPDATE MARKERS LIST
    {
        markers.list <- list()
        selected.markers <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_2")
        {
            if(length(input[["t_2_file_list"]])>0 && input[["t_2_file_list"]] != "")
            {
                f.name <- input[["t_2_file_list"]]
                if(length(app.variables$fcs.files[[f.name]])>0)
                {
                    if(app.variables$fcs.files[[f.name]][["set"]] ==  input[["t_2_set_list"]])
                    {
                        markers.list <- app.variables$fcs.files[[f.name]][["markers"]]
                    }
                }
            }
            selected.markers <- markers.list
            if(!is.na(input$t_2_markers_list) && length(input$t_2_markers_list)>0)
            {
                selected.markers <- as.list(unlist(input$t_2_markers_list)[input$t_2_markers_list%in%markers.list])
            }
            if(length(selected.markers)>0)
            {
                if(length(app.variables$fcs.files[[f.name]])>0)
                {
                    pop.col <- app.variables$fcs.files[[f.name]][["populations_column"]]
                    pop.col <- app.variables$fcs.files[[f.name]][["markers"]][[pop.col]]
                    #==
                    selected.markers <- as.list(unlist(selected.markers)[which(selected.markers!=pop.col)])
                }
            }
        }
        updateSelectInput(session, "t_2_markers_list", "Select Markers", choices = markers.list, selected = selected.markers)
    })
    
    observe( #UPDATE VISUALIZED MARKERS LIST
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_2")
        {
            markers.list <- list()
            selected.markers.1 <- NULL
            selected.markers.2 <- NULL
            if(length(input[["t_2_file_list"]])>0 && input[["t_2_file_list"]] != "")
            {
                if(length(input[["t_2_markers_list"]])>0 && input[["t_2_markers_list"]]!="")
                {
                    f.name <- input[["t_2_file_list"]]
                    if(length(app.variables$fcs.files[[f.name]])>0)
                    {
                        markers.list <- input[["t_2_markers_list"]]
                    }
                }
            }
            if(!is.na(input$t_2_m1) && input$t_2_m1!="")
            {
                selected.markers.1 <- as.list(unlist(input$t_2_m1)[input$t_2_m1%in%markers.list])
            }
            else
            {
                selected.markers.1 <- markers.list
            }
            if(!is.na(input$t_2_m2) && input$t_2_m2!="")
            {
                selected.markers.2 <- as.list(unlist(input$t_2_m2)[input$t_2_m2%in%markers.list])
            }
            else
            {
                selected.markers.2 <- markers.list
            }
            
            updateSelectInput(session, "t_2_m1", "Select 1st Marker", choices = markers.list, selected = selected.markers.1)
            updateSelectInput(session, "t_2_m2", "Select 2nd Marker", choices = markers.list, selected = selected.markers.2)
        }
    })
    
    output$t_2_means <- renderUI(
    {
        rendered.UI <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_2")
        { 
            if(length(input[["t_2_file_list"]])>0 && input[["t_2_file_list"]] != "")
            {
                if(length(input[["t_2_markers_list"]])>0)
                {
                    if(length(input[["t_2_pop_list"]])>0 && input[["t_2_pop_list"]]!="")
                    {
                        f.name <- input[["t_2_file_list"]]
                        if(length(app.variables$fcs.files[[f.name]])>0)
                        {
                            tmp.mat <- app.variables$output.matrices[[f.name]]
                            markers.list <- input[["t_2_markers_list"]]
                            current.pop <- input$t_2_pop_list
                            
                            min.value <- as.numeric(input[["t_2_min_value"]])
                            max.value <- as.numeric(input[["t_2_max_value"]])
                            
                            rendered.UI <- lapply(1:length(markers.list), function(j)
                            {
                                tmp.UI <- NULL
                                used.value <- min.value
                                if(length(app.variables$fcs.files[[f.name]][["markers_default_values"]][[current.pop]])>0)
                                {
                                    used.value <- app.variables$fcs.files[[f.name]][["markers_default_values"]][[current.pop]][[1]][[j]]
                                }
                                if(length(input[[paste0("t_2_mean_",j)]])>0 && input[[paste0("t_2_mean_",j)]]!="")
                                {
                                    used.value <- input[[paste0("t_2_mean_",j)]]
                                }
                                marker.name <- markers.list[[j]]
                                marker.id <- which(app.variables$fcs.files[[f.name]][["markers"]] == marker.name)
                                if(length(marker.id)>0)
                                {
                                    marker.id <- as.numeric(marker.id[[1]])
                                    mean.step <- trunc(1000*((min.value+max.value)/100))/1000
                                    tmp.UI <- sliderInput(paste0("t_2_mean_",j), marker.name,
                                                          min = min.value,
                                                          max = max.value,
                                                          step=mean.step,
                                                          value=used.value,
                                                          width = "80%")
                                }
                                return(tmp.UI)
                            })
                        }
                    }
                }
            }
        }
    })
    
    output$t_2_sd <- renderUI(
    {
        rendered.UI <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_2")
        { 
            if(length(input[["t_2_file_list"]])>0 && input[["t_2_file_list"]] != "")
            {
                if(length(input[["t_2_markers_list"]])>0)
                {
                    if(length(input[["t_2_pop_list"]])>0 && input[["t_2_pop_list"]]!="")
                    {
                        f.name <- input[["t_2_file_list"]]
                        if(length(app.variables$fcs.files[[f.name]])>0)
                        {
                            tmp.mat <- app.variables$output.matrices[[f.name]]
                            markers.list <- input[["t_2_markers_list"]]
                            current.pop <- input$t_2_pop_list
                            
                            min.value <- as.numeric(input[["t_2_min_value"]])
                            max.value <- as.numeric(input[["t_2_max_value"]])
                            
                            rendered.UI <- lapply(1:length(markers.list), function(j)
                            {
                                tmp.UI <- NULL
                                used.value <- min.value
                                if(length(app.variables$fcs.files[[f.name]][["markers_default_values"]][[current.pop]])>0)
                                {
                                    used.value <- app.variables$fcs.files[[f.name]][["markers_default_values"]][[current.pop]][[2]][[j]]
                                }
                                if(length(input[[paste0("t_2_sd_",j)]])>0 && input[[paste0("t_2_sd_",j)]]!="")
                                {
                                    used.value <- input[[paste0("t_2_sd_",j)]]
                                }
                                marker.name <- markers.list[[j]]
                                marker.id <- which(app.variables$fcs.files[[f.name]][["markers"]] == marker.name)
                                if(length(marker.id)>0)
                                {
                                    marker.id <- as.numeric(marker.id[[1]])
                                    sd.max <- (min.value+max.value)/5
                                    sd.step <- trunc(1000*sd.max/100)/1000
                                    tmp.UI <- sliderInput(paste0("t_2_sd_",j), marker.name,
                                                          min = 0,
                                                          max = sd.max,
                                                          step=sd.step,
                                                          value=used.value,
                                                          width = "80%")
                                }
                                return(tmp.UI)
                            })
                        }
                    }
                }
            }
        }
    })
    
    observe( #UPDATE POPULATIONS LIST
    {
        pop.list <- list()
        selected.pop <- list()
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_2")
        {
            if(length(input[["t_2_set_list"]])>0)
            {
                if(length(input[["t_2_file_list"]])>0 && input[["t_2_file_list"]] != "")
                {
                    f.name <- input[["t_2_file_list"]]
                    if(length(app.variables$fcs.files[[f.name]])>0 && app.variables$fcs.files[[f.name]][["set"]] ==  input[["t_2_set_list"]])
                    {
                        f.id <- which(names(app.variables$fcs.files)==f.name)[[1]]
                        fcs.temp <- app.variables$fcs.files[[f.name]]
                        fcs.mat <- app.variables$output.matrices[[f.name]]
                        pop.list <- list()
                        for(i in 1:length(app.variables$populations.list[[f.name]]))
                        {
                            if(length(app.variables$populations.list[[f.name]][[i]])>0)
                            {
                                pop.list <- list(unlist(pop.list), names(app.variables$populations.list[[f.name]])[i])
                            }
                        }
                    }
                }
                pop.list <- unlist(pop.list)
                selected.pop <- pop.list
                if(!is.na(input$t_2_pop_list) && input$t_2_pop_list!="")
                {
                    selected.pop <- as.list(unlist(input$t_2_pop_list)[input$t_2_pop_list%in%pop.list])
                }
            }
        }
        updateSelectInput(session, "t_2_pop_list", "Select a Population", choices=pop.list, selected = selected.pop)
    })
    
    observe( #UPDATE MARKERS DEFAULT VALUES
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_2")
        {
            if(length(input[["t_2_file_list"]])>0 && input[["t_2_file_list"]] != "")
            {
                f.name <- input[["t_2_file_list"]]
                if(length(app.variables$fcs.files[[f.name]])>0)
                {
                    tmp.mat <- app.variables$output.matrices[[f.name]]
                    if(length(app.variables$populations.list[[f.name]])>0)
                    {
                        lapply(1:length(app.variables$populations.list[[f.name]]), function(i)
                        {
                            if(length(app.variables$populations.list[[f.name]][[i]])>0)
                            {
                                pop <- tmp.mat[as.numeric(app.variables$populations.list[[f.name]][[i]]),]
                                if(length(app.variables$populations.list[[f.name]][[i]])==1)
                                {
                                    pop <- t(tmp.mat[as.numeric(app.variables$populations.list[[f.name]][[i]]),])
                                }
                                pop.name <- names(app.variables$populations.list[[f.name]])[i]
                                app.variables$fcs.files[[f.name]][["markers_default_values"]][[pop.name]] <-  extract.position.from.pop(pop)
                            }
                        })
                    }
                }
            }
        }
    })
    
    observe( #UPDATE SD ET MEAN SLIDERS VALUES
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_2")
        {
            if(length(input[["t_2_file_list"]])>0 && input[["t_2_file_list"]] != "")
            {
                if(length(input[["t_2_pop_list"]])>0 && input[["t_2_pop_list"]] != "")
                {
                    f.name <- input[["t_2_file_list"]]
                    if(length(app.variables$fcs.files[[f.name]])>0 && 
                       length(app.variables$fcs.files[[f.name]][["markers_default_values"]])>0)
                    {
                        markers.values <- app.variables$fcs.files[[f.name]][["markers_default_values"]][[input[["t_2_pop_list"]]]]
                        lapply(1:length(markers.values[[1]]), function(i)
                        {
                            updateSliderInput(session, paste0("t_2_mean_",i), value = markers.values[[1]][[i]])
                            updateSliderInput(session, paste0("t_2_sd_",i), value = markers.values[[2]][[i]])
                        })
                    }
                }
            }
        }
    })
    
    observeEvent(input$t_2_pop_new, #ADD NEW POP
    {
        shinyjs::disable("t_2_pop_new")
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_2")
        {
            if(length(input[["t_2_file_list"]])>0 && input[["t_2_file_list"]] != "" && 
               app.variables$fcs.files[[input[["t_2_file_list"]]]][["set"]] ==  input[["t_2_set_list"]])
            {
                if(length(input[["t_2_pop_events"]])>0 && input[["t_2_pop_events"]] != "")
                {
                    progress <- Progress$new()
                    progress$set(message = "Adding Population", value = 0)
                    
                    fcs.id <- which(names(app.variables$fcs.files) == input[["t_2_file_list"]])[[1]]
                    
                    tmp.matrix <- matrix(nrow=as.numeric(input[["t_2_pop_events"]]), 
                                         ncol=ncol(app.variables$output.matrices[[input[["t_2_file_list"]]]]))
                    
                    pop.col <- app.variables$fcs.files[[fcs.id]][["populations_column"]]
                    for(i in (1:ncol(tmp.matrix))[-pop.col])
                    {
                        m.orig <- runif(1,
                                        as.numeric(input$t_2_min_value),
                                        as.numeric(input$t_2_max_value))
                        sd.orig <- runif(1,0,1)*1.76
                        
                        tmp.matrix[,i] <- rtruncnorm(as.numeric(input[["t_2_pop_events"]]),
                                                     as.numeric(input$t_2_min_value),
                                                     as.numeric(input$t_2_max_value),
                                                     m.orig,
                                                     sd.orig)
                        
                        progress$inc(1/(ncol(tmp.matrix)-1), detail = "Creating points")
                    }
                    tmp.matrix[,pop.col] <- max(app.variables$output.matrices[[fcs.id]][,pop.col])+1
                    first.id <- nrow(app.variables$output.matrices[[fcs.id]])+1
                    
                    app.variables$output.matrices[[fcs.id]] <<-
                        rbind(app.variables$output.matrices[[fcs.id]],
                              tmp.matrix)
                    app.variables$populations.list[[fcs.id]][[length(app.variables$populations.list[[fcs.id]])+1]] <<- 
                        first.id:(first.id+nrow(tmp.matrix)-1)
                    names(app.variables$populations.list[[fcs.id]])[[length(app.variables$populations.list[[fcs.id]])]] <<- 
                        length(app.variables$populations.list[[fcs.id]])
                    app.variables$used.events[[fcs.id]] <<- c(unlist(app.variables$used.events[[fcs.id]]), rep(T,nrow(tmp.matrix)))
                    
                    progress$set(message = "Population added", value = 1)
                    delay(700, progress$close())
                }
            }
        }
        delay(500, shinyjs::enable("t_2_pop_new"))
    })
    
    observeEvent(input$t_2_pop_inc, #INCREASE POP SIZE
    {
        shinyjs::disable("t_2_pop_inc")
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_2")
        {
            if(length(input[["t_2_file_list"]])>0 && input[["t_2_file_list"]] != "" && 
               length(app.variables$fcs.files[[input[["t_2_file_list"]]]])>0 && 
               app.variables$fcs.files[[input[["t_2_file_list"]]]][["set"]] ==  input[["t_2_set_list"]])
            {
                if(length(input[["t_2_pop_events"]])>0 && input[["t_2_pop_events"]] != "")
                {
                    if(length(input[["t_2_pop_list"]])>0 && input[["t_2_pop_list"]] != "")
                    {
                        
                        f.name <- input[["t_2_file_list"]]
                        pop.list <- app.variables$populations.list[[f.name]]
                        current.pop <- as.numeric(which(names(pop.list)%in%input[["t_2_pop_list"]])[[1]])
                        
                        progress <- Progress$new()
                        progress$set(message = paste0("Increasing size of ", input[["t_2_pop_list"]]), value = 0)
                        
                        #AUGMENTATION DE LA TAILLE DE LA POP
                        pop.col <- app.variables$fcs.files[[f.name]][["populations_column"]]
                        if(length(pop.list[[current.pop]])>0)
                        {
                            pop <- app.variables$output.matrices[[f.name]][ pop.list[[current.pop]], ]
                            inc.coef <- 100*as.numeric(input[["t_2_pop_events"]])/nrow(pop)
                            new.mat <- create.pop.from.pop(pop, inc.coef, unused.columns = pop.col, limited=T)
                            first.id <- nrow(app.variables$output.matrices[[f.name]])
                            
                            app.variables$output.matrices[[f.name]] <<- rbind(app.variables$output.matrices[[f.name]],
                                                                              new.mat)
                            app.variables$populations.list[[f.name]][[current.pop]] <<-
                                c(unlist(pop.list[[current.pop]]), (first.id+1):(first.id+nrow(new.mat)))
                            app.variables$used.events[[f.name]] <<- c(app.variables$used.events[[f.name]],
                                                                      rep(T,nrow(new.mat)))
                        }
                        
                        progress$inc(1, detail = "done")
                        progress$set(message = "Size increased", value = 1)
                        
                        delay(700, progress$close())
                    }
                }
            }
        }
        delay(500, shinyjs::enable("t_2_pop_inc"))
    })
    
    observeEvent(input$t_2_pop_red, #DECREASE POP SIZE
    {
        shinyjs::disable("t_2_pop_red")
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_2")
        {
            if(length(input[["t_2_file_list"]])>0 && input[["t_2_file_list"]] != "" && 
               length(app.variables$fcs.files[[input[["t_2_file_list"]]]])>0 && 
               app.variables$fcs.files[[input[["t_2_file_list"]]]][["set"]] ==  input[["t_2_set_list"]])
            {
                if(length(input[["t_2_pop_events"]])>0 && input[["t_2_pop_events"]] != "")
                {
                    if(length(input[["t_2_pop_list"]])>0 && input[["t_2_pop_list"]] != "")
                    {
                        f.name <- input[["t_2_file_list"]]
                        pop.list <- app.variables$populations.list[[f.name]]
                        current.pop <- as.numeric(which(names(pop.list)%in%input[["t_2_pop_list"]])[[1]])
                        
                        progress <- Progress$new()
                        progress$set(message = paste0("Reducing size of ", input[["t_2_pop_list"]]), value = 0)
                        
                        #AUGMENTATION DE LA TAILLE DE LA POP
                        pop.col <- app.variables$fcs.files[[f.name]][["populations_column"]]
                        
                        if(length(pop.list[[current.pop]])>0)
                        {
                            pop <- app.variables$output.matrices[[f.name]][ pop.list[[current.pop]], ]
                            red.coef <- 100*as.numeric(input[["t_2_pop_events"]])/nrow(pop)
                            if(red.coef>100)
                            {
                                red.coef <- 100
                            }
                            removed.events <- reduce.population(pop.list[[current.pop]], red.coef)
                            removed.events.ids <- which(pop.list[[current.pop]]%in%removed.events)
                            
                            app.variables$used.events[[f.name]][removed.events] <<- F
                            app.variables$populations.list[[f.name]][[current.pop]] <<- pop.list[[current.pop]][-removed.events.ids]
                        }
                        
                        progress$inc(1, detail = "done")
                        progress$set(message = "Size decreased", value = 1)
                        delay(700, progress$close())
                    }
                }
            }
        }
        delay(500, shinyjs::enable("t_2_pop_red"))
    })
    
    output$t_2_plot <- renderPlot( #RENDER PLOT
    {
        output.plot <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_2")
        {
            if(length(input[["t_2_file_list"]])>0 && input[["t_2_file_list"]] != "")
            {
                if(length(input[["t_2_m1"]])>0 && input[["t_2_m1"]] != "")
                {
                    if(length(input[["t_2_m2"]])>0 && input[["t_2_m2"]] != "")
                    {
                        if(length(input[["t_2_pop_list"]])>0 && input[["t_2_pop_list"]] != "")
                        {
                            f.name <- isolate(input[["t_2_file_list"]])
                            if(length(app.variables$fcs.files[[f.name]])>0 &&
                               app.variables$fcs.files[[f.name]][["set"]] ==  input[["t_2_set_list"]])
                            {
                                used.events <- app.variables$used.events[[f.name]]
                                fcs.mat <- app.variables$output.matrices[[f.name]][used.events, ]
                                f.id <- which(names(app.variables$fcs.files)==f.name)[[1]]
                                
                                m1.id <- which(app.variables$fcs.files[[f.name]][["markers"]]==input[["t_2_m1"]])[[1]]
                                m2.id <- which(app.variables$fcs.files[[f.name]][["markers"]]==input[["t_2_m2"]])[[1]]
                                
                                pop.num <- rep("Other",nrow(app.variables$output.matrices[[f.name]]))
                                pop.events <- as.numeric(unlist(app.variables$populations.list[[f.name]][[input[["t_2_pop_list"]]]]))
                                if(length(pop.events)>0)
                                {
                                    pop.num[pop.events] <- "Selected Population"
                                }
                                pop.num <- pop.num[used.events]
                                tmp.dataframe <- data.frame(P1=unlist(fcs.mat[,as.integer(m1.id)]),
                                                            P2=unlist(fcs.mat[,as.integer(m2.id)]),
                                                            pop=unlist(pop.num))
                                tmp.dataframe$pop <- as.factor(tmp.dataframe$pop)
                                
                                output.plot <- ggplot(tmp.dataframe, aes(x=P1, y=P2, color=pop)) + 
                                    geom_point(size = 0.1, stroke = 0, shape = 16) +
                                    xlab(input[["t_2_m1"]]) + 
                                    ylab(input[["t_2_m2"]]) +
                                    theme(legend.position = "top") +
                                    xlim(as.numeric(input[["t_2_min_value"]]), as.numeric(input[["t_2_max_value"]])) +
                                    ylim(as.numeric(input[["t_2_min_value"]]), as.numeric(input[["t_2_max_value"]])) 
                            }
                        }
                    }
                }
            }
        }
        return(output.plot)
    })
    
    observeEvent(input$t_2_pop_move, #MOVE POPULATION AVEC SLIDERS
    {
        shinyjs::disable("t_2_pop_move")
        
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_2")
        {
            if(length(input[["t_2_file_list"]])>0 && input[["t_2_file_list"]] != "")
            {
                if(length(input[["t_2_pop_list"]])>0 && input[["t_2_pop_list"]] != "")
                {
                    progress <- Progress$new()
                    progress$set("Moving Population", value=0)
                    
                    f.name <- input[["t_2_file_list"]]
                    if(length(app.variables$fcs.files[[f.name]])>0 && app.variables$fcs.files[[f.name]][["set"]] ==  input[["t_2_set_list"]])
                    {
                        f.id <- which(names(app.variables$fcs.files)==f.name)[[1]]
                        
                        pop.column <- app.variables$fcs.files[[f.name]][["populations_column"]]
                        markers.values <- app.variables$fcs.files[[f.name]][["markers_default_values"]][[input[["t_2_pop_list"]]]]
                        markers.list <- (1:length(markers.values[[1]]))[-pop.column]
                        events.list <- app.variables$populations.list[[f.name]][[input[["t_2_pop_list"]]]]
                        if(length(events.list)>0)
                        {
                            current.pop <- app.variables$output.matrices[[f.name]][events.list,-pop.column]
                            mean.list <- c()
                            sd.list <- c()
                            
                            for(i in 1:length(markers.values[[1]]))
                            {
                                if(length(input[[paste0("t_2_mean_",i)]])>0 && input[[paste0("t_2_mean_",i)]] != "")
                                {
                                    mean.list[[i]] <- as.numeric(input[[paste0("t_2_mean_",i)]])
                                    sd.list[[i]] <- as.numeric(input[[paste0("t_2_sd_",i)]])
                                }
                                else
                                {
                                    mean.list[[i]] <- markers.values[[1]][[i]]
                                    sd.list[[i]] <- markers.values[[2]][[i]]
                                }
                                progress$inc(1/length(markers.values[[1]]), 
                                             detail=paste0(app.variables$fcs.files[[f.name]][["markers"]][[i]], " changed"))
                            }
                            if(length(mean.list)>0)
                            {
                                app.variables$output.matrices[[f.name]][events.list,markers.list] <<- 
                                    move.population(current.pop, mean.list, sd.list, limited=T, 
                                                    max.val = as.numeric(input[["t_2_max_value"]]), 
                                                    min.val = as.numeric(input[["t_2_min_value"]]))
                            }
                        }
                    }
                    progress$set(paste0(input[["t_2_pop_list"]], " moved"), value=1)
                    delay(500, progress$close())
                }
            }
        }
        
        
        delay(500, shinyjs::enable("t_2_pop_move"))
    })
    
    observeEvent(input$t_2_pop_remove, #DEL POPULATION
    {
        shinyjs::disable("t_2_pop_remove")
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_2")
        {
            if(length(input[["t_2_file_list"]])>0 && input[["t_2_file_list"]] != "")
            {
                if(length(input[["t_2_pop_list"]])>0 && input[["t_2_pop_list"]] != "")
                {
                    progress <- Progress$new()
                    progress$set("Deleting Population", value=0)
                    f.name <- input[["t_2_file_list"]]
                    if(length(app.variables$fcs.files[[f.name]])>0 && app.variables$fcs.files[[f.name]][["set"]] ==  input[["t_2_set_list"]])
                    {
                        f.id <- which(names(app.variables$fcs.files)==f.name)[[1]]
                        
                        pop.column <- app.variables$fcs.files[[f.name]][["populations_column"]]
                        
                        events.list <- app.variables$populations.list[[f.name]][[input[["t_2_pop_list"]]]]
                        if(length(events.list)>0)
                        {
                            app.variables$used.events[[f.name]][events.list] <<- F
                            app.variables$populations.list[[f.name]][[input[["t_2_pop_list"]]]] <<- list()
                        }
                    }
                    progress$set(paste0(input[["t_2_pop_list"]], " removed"), value=1)
                    delay(500, progress$close())
                }
            }
        }
        delay(500, shinyjs::enable("t_2_pop_remove"))
    })
    
    
    
    
    
    
    
    
    
    #========================================GENERATE CONTROL FILES================================================
    
    
    
    
    
    
    
    
    observe( #CLEAR UI
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_3")
        {
            if(length(input[["t_3_file_list"]])>0 && input[["t_3_file_list"]] != "")
            {
                f.name <- input[["t_3_file_list"]]
                if(length(app.variables$fcs.files[[f.name]])>0)
                {
                    shinyjs::show("t_3_2")
                    shinyjs::show("t_3_3")
                    shinyjs::enable("t_3_control_generate")
                }
                else
                {
                    hide("t_3_2")
                    hide("t_3_3")
                    shinyjs::disable("t_3_control_generate")
                }
            }
            else
            {
                hide("t_3_2")
                hide("t_3_3")
                shinyjs::disable("t_3_control_generate")
            }
        }
        else
        {
            hide("t_3_2")
            hide("t_3_3")
            shinyjs::disable("t_3_control_generate")
        }
    })
    
    observe( #UPDATE SETS LIST
    {
        update.sets.list(3)
    })
    
    observe( #UPDATE FILES LIST
    {
        update.files.list(3)
    })
    
    observe( #UPDATE MARKERS LIST
    {
        markers.list <- list()
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(length(input[["t_3_file_list"]])>0 && input[["t_3_file_list"]] != "")
            {
                f.name <- input[["t_3_file_list"]]
                if(length(app.variables$fcs.files[[f.name]]) > 0 && app.variables$fcs.files[[f.name]][["set"]] ==  input[["t_3_set_list"]])
                {
                    f.id <- which(names(app.variables$fcs.files)==f.name)[[1]]
                    markers.list <- unlist(app.variables$fcs.files[[f.name]][["markers"]])
                    
                    pop.col <- app.variables$fcs.files[[f.name]][["populations_column"]]
                    markers.list <- markers.list[-pop.col]
                }
            }
        }
        updateSelectInput(session, "t_3_var_markers_list", "Select 1st Marker", choices = markers.list, selected = markers.list)
    })
    
    observe( #UPDATE POPULATIONS LIST
    {
        pop.list <- list()
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(length(input[["t_3_set_list"]])>0)
            {
                if(length(input[["t_3_file_list"]])>0 && input[["t_3_file_list"]] != "")
                {
                    f.name <- input[["t_3_file_list"]]
                    
                    if(length(app.variables$fcs.files[[f.name]]) > 0 && app.variables$fcs.files[[f.name]][["set"]] ==  input[["t_3_set_list"]] &&
                       length(app.variables$populations.list[[f.name]])>0)
                    {
                        fcs.temp <- app.variables$fcs.files[[f.name]]
                        fcs.mat <- app.variables$output.matrices[[f.name]]
                        
                        pop.list <- list()
                        for(i in 1:length(app.variables$populations.list[[f.name]]))
                        {
                            if(length(app.variables$populations.list[[f.name]][[i]])>0)
                            {
                                pop.list <- list(unlist(pop.list), names(app.variables$populations.list[[f.name]])[i])
                            }
                        }
                        pop.list <- unlist(pop.list)
                    }
                }
            }
        }
        updateSelectInput(session, "t_3_var_pop_list", "Select a Population", choices=pop.list, selected = pop.list)
    })
    
    observeEvent(input$t_3_control_generate,#GENERATE CONTROL FILES
    {
        shinyjs::disable("t_3_control_generate")
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(length(input[["t_3_set_list"]])>0)
            {
                if(length(input[["t_3_file_list"]])>0 && input[["t_3_file_list"]] != "")
                {
                    if(length(input[["t_3_var_markers_list"]])>0 && input[["t_3_var_markers_list"]] != "")
                    {
                        if(length(input[["t_3_var_pop_list"]])>0 && input[["t_3_var_pop_list"]] != "")
                        {
                            f.name <- input[["t_3_file_list"]]
                            if(length(app.variables$fcs.files[[f.name]])> 0 && app.variables$fcs.files[[f.name]][["set"]] ==  input[["t_3_set_list"]])
                            {
                                progress <- Progress$new()
                                progress$set(message = "Generating Control Files", value = 0)
                                update.log("GENERATING CONTROL FILES")
                                
                                used.events <- app.variables$used.events[[f.name]]
                                f.set <- app.variables$fcs.files[[f.name]][["set"]]
                                f.mat <- app.variables$output.matrices[[f.name]]
                                f.populations <- app.variables$populations.list[[f.name]]
                                update.log(paste("==== FROM:", f.name))
                                
                                nmb.files <- as.numeric(input[["t_3_nmb_ctrl"]])
                                unused.markers <- as.numeric(unlist(which(!(app.variables$fcs.files[[f.name]][["markers"]]%in%input[["t_3_var_markers_list"]]))))
                                selected.populations <- as.numeric(unlist(which(names(app.variables$populations.list[[f.name]])%in%input[["t_3_var_pop_list"]])))
                                
                                
                                update.log(paste("==== MODIFIED POPULATIONS:",  list.to.string(input[["t_3_var_pop_list"]])))
                                update.log(paste("==== VARIABLE MARKERS:",  list.to.string(input[["t_3_var_markers_list"]])))
                                
                                
                                lapply(1:nmb.files, function(curr.file)
                                {
                                    tmp.mat <- f.mat
                                    lapply(selected.populations, function(curr.pop)
                                    {
                                        if(length(f.populations[[curr.pop]])>0)
                                        {
                                            if(length(f.populations[[curr.pop]])==1)
                                            {
                                                pop <- t(tmp.mat[f.populations[[curr.pop]],])
                                                pop.position <- extract.position.from.pop(pop, unused.markers)
                                                
                                                tmp.mat[f.populations[[curr.pop]], ] <<-
                                                    t(move.population(pop, pop.position[[1]], pop.position[[2]],
                                                                      limited = T, unused.columns = unused.markers))
                                            }
                                            else
                                            {
                                                
                                                pop <- tmp.mat[f.populations[[curr.pop]],]
                                                pop.position <- extract.position.from.pop(pop, unused.markers)
                                                
                                                tmp.mat[f.populations[[curr.pop]], ] <<-
                                                    move.population(pop, pop.position[[1]], pop.position[[2]],
                                                                    limited = T, unused.columns = unused.markers)
                                            }
                                            
                                        }
                                    })
                                    
                                    tmp.name <<- paste0("CTRL_",f.name,"_", curr.file, "_", trunc(as.numeric(Sys.time())))
                                    app.variables$output.matrices[[tmp.name]] <<- as.matrix(tmp.mat)
                                    app.variables$used.events[[tmp.name]] <<- used.events
                                    app.variables$populations.list[[tmp.name]] <<- isolate(app.variables$populations.list[[f.name]])
                                    app.variables$fcs.files[[tmp.name]] <<- isolate(app.variables$fcs.files[[f.name]])
                                    app.variables$fcs.files[[tmp.name]][["file"]]@exprs <<- as.matrix(tmp.mat)
                                    app.variables$fcs.files[[tmp.name]][["name"]] <<- tmp.name
                                    app.variables$fcs.files[[tmp.name]][["type"]] <<- isolate(app.variables$fcs.files[[f.name]][["type"]])
                                    app.variables$fcs.files[[tmp.name]][["source_ctrl"]] <<- input[["t_3_set_list"]]
                                    progress$inc(1/nmb.files, detail=paste0("CTRL ", curr.file," generated"))
                                    update.log(paste("========", tmp.name, "generated"))
                                })
                                
                                
                                progress$set(message = "Control Files Generated", value = 1)
                                update.log("CONTROL FILES GENERATED")
                                update.log("")
                                delay(500, progress$close())
                            }
                        }
                    }
                }
            }
        }
        delay(500, shinyjs::enable("t_3_control_generate"))
    })
    
    
    
    
    
    
    
    #===========================================GENERATE MUTANTS===================================================
    
    
    
    
    
    
    observe( #CLEAR UI
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_6")
        {
            if(length(input[["t_6_files_list"]])>0 && input[["t_6_files_list"]] != "")
            {
                f.name <- input[["t_6_files_list"]]
                if(length(app.variables$fcs.files[[f.name]])>0)
                {
                    shinyjs::enable("t_6_mut_generate")
                }
                else
                {
                    shinyjs::disable("t_6_mut_generate")
                }
            }
            else
            {
                shinyjs::disable("t_6_mut_generate")
            }
        }
        else
        {
            shinyjs::disable("t_6_mut_generate")
        }
    })
    
    observe( #UPDATE SETS LIST
    {
        update.sets.list(6)
    })

    observe( #UPDATE FILES LIST
    {
        files.list <- list()
        selected.file <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(!is.na(input[["t_6_set_list"]])  && input[["t_6_set_list"]]!="")
            {
                for(i in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[i]])>0 && app.variables$fcs.files[[i]][["set"]] ==  input[["t_6_set_list"]] &&
                       app.variables$fcs.files[[i]][["type"]] == "CTRL")
                    {
                        files.list[[length(files.list)+1]] <- names(app.variables$fcs.files)[[i]]
                    }
                }
            }
            selected.file <- files.list
            if(!is.na(input$t_6_files_list) && input$t_6_files_list!="")
            {
                selected.file <- input$t_6_files_list
            }
        }
        updateSelectInput(session, "t_6_files_list", "Select Control Files", choices=files.list, selected=selected.file)
    })
    
    output$t_6_1 <- renderUI( #RENDER UI POUR LE CTRL
    {
        rendered.ui <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(length(input$t_6_files_list)>0 && input$t_6_files_list != "")
            {
                current.ctrl <- input$t_6_files_list
                #==
                if(length(app.variables$fcs.files[[current.ctrl]])>0 &&
                   app.variables$fcs.files[[current.ctrl]][["set"]] ==  input[["t_6_set_list"]])
                {
                    mutant.list <- list()
                    for(i in 1:length(app.variables$fcs.files))
                    {
                        if(length(app.variables$fcs.files[[i]])>0 && app.variables$fcs.files[[i]][["set"]] ==  input[["t_6_set_list"]])
                        {
                            if(app.variables$fcs.files[[i]][["set"]] ==  input[["t_6_set_list"]] && app.variables$fcs.files[[i]][["type"]] == "MUT" &&
                               app.variables$fcs.files[[i]][["source_ctrl"]] == current.ctrl)
                            {
                                mutant.list[[length(mutant.list)+1]] <- names(app.variables$fcs.files)[[i]]
                            }
                        }
                    }
                    #=
                    pop.list <- app.variables$populations.list[[current.ctrl]]
                    pop.tags <- list()
                    if(length(mutant.list)>0 && !is.null(app.variables$loaded.mutant))
                    {
                        if(length(app.variables$reduction.percentages[[app.variables$loaded.mutant]])>0)
                        {
                            pop.tags <- lapply(1:length(pop.list), function(current.pop)
                            {
                                return(tagList(
                                    column(
                                        width=4,
                                        h4(names(pop.list)[[current.pop]])
                                    ),
                                    column(
                                        width=8,
                                        numericInput(paste0("t_6_val_",current.pop), "Size Reduction %",
                                                     value=app.variables$reduction.percentages[[app.variables$loaded.mutant]][[current.pop]])
                                    )
                                ))
                            })
                        }
                        else
                        {
                            pop.tags <- lapply(1:length(pop.list), function(current.pop)
                            {
                                return(tagList(
                                    column(
                                        width=4,
                                        h4(current.pop)
                                    ),
                                    column(
                                        width=8,
                                        numericInput("t_6_val_", "Reduction %", value=0)
                                    )
                                ))
                            })
                        }
                        
                    }
                    #=
                    rendered.ui <- box(
                        title=app.variables$fcs.files[[input$t_6_files_list]][["name"]],
                        width=12, collapsible=T,
                        column(
                            width=7,
                            div(
                                width="90%",style="padding:0",
                                column(
                                    width=9,style="padding:0",
                                    selectInput("t_6_mutants_list", NULL, choices = mutant.list, width="100%")
                                ),
                                column(
                                    width=2,style="padding:0;margin-left:1%",
                                    actionButton("t_6_load_mutant", "Load", width="100%")
                                )
                            )
                        ),
                        column(
                            width=5,
                            pop.tags
                        )
                    )
                }
            }
        }
        return(rendered.ui)
    })
    
    observeEvent(input[["t_6_load_mutant"]], #OBSERVE EVENT DU BOUTON LOAD MUTANT
    {
        shinyjs::disable("t_6_load_mutant")
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(length(input$t_6_files_list)>0 && input$t_6_files_list != "")
            {
                if(length(app.variables$fcs.files[[input$t_6_files_list]])>0 &&
                   app.variables$fcs.files[[input$t_6_files_list]][["set"]] ==  input[["t_6_set_list"]])
                {
                    progress <- Progress$new()
                    progress$set(message = "Loading Mutant", value = 0)
                    current.mut <- input[["t_6_mutants_list"]]
                    app.variables$loaded.mutant <<- current.mut
                    
                    pop.list <- app.variables$populations.list[[input$t_6_files_list]]
                    if(length(pop.list)>0)
                    {
                        if(length(app.variables$reduction.percentages[[current.mut]])>0)
                        {
                            sapply(1:length(pop.list), function(j)
                            {
                                updateNumericInput(session, paste0("t_6_val_",j), "Reduction %",
                                                   value=app.variables$reduction.percentages[[current.mut]][[j]])
                            })
                        }
                    }
                    
                    progress$set(message = "Mutant loaded", value = 1)
                    delay(500, progress$close())
                }
            }
        }
        delay(500, shinyjs::enable("t_6_load_mutant"))
    })
    
    add.mutant <- function(ctrl.name)
    {
        shinyjs::disable("t_6_mut_add")
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(length(ctrl.name)>0 && ctrl.name != "")
            {
                if(length(app.variables$fcs.files[[ctrl.name]]) > 0 &&
                   app.variables$fcs.files[[ctrl.name]][["set"]] ==  input[["t_6_set_list"]])
                {
                    progress <- Progress$new()
                    progress$set(message = "Adding Mutant", value = 0)
                    
                    f.name <- app.variables$fcs.files[[ctrl.name]][["name"]]
                    tmp.name <- paste0("MUT_",f.name,"_",trunc(as.numeric(Sys.time())))
                    app.variables$output.matrices[[tmp.name]] <<- isolate(app.variables$output.matrices[[ctrl.name]])
                    app.variables$populations.list[[tmp.name]] <<- isolate(app.variables$populations.list[[ctrl.name]])
                    app.variables$fcs.files[[tmp.name]] <<- isolate(app.variables$fcs.files[[ctrl.name]])
                    app.variables$fcs.files[[tmp.name]][["type"]] <<- "MUT"
                    app.variables$fcs.files[[tmp.name]][["source_ctrl"]] <<- ctrl.name
                    app.variables$fcs.files[[tmp.name]][["name"]] <<- tmp.name
                    app.variables$fcs.files[[tmp.name]][["file"]]@exprs <<- isolate(app.variables$output.matrices[[ctrl.name]])
                    app.variables$reduction.percentages[[tmp.name]] <<- rep(0,length(app.variables$populations.list[[ctrl.name]]))
                    app.variables$used.events[[tmp.name]] <<- isolate(app.variables$used.events[[ctrl.name]])
                    
                    progress$set(message = "Mutant added", value = 1)
                    update.log(paste(f.name, input$t_6_set_list))
                    
                    delay(500, progress$close())
                }
            }
        }
        delay(500, shinyjs::enable("t_6_mut_add"))
    }
    
    observeEvent(input[["t_6_mut_add"]], #OBSERVE EVENT DU BOUTON ADD MUTANT
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(length(app.variables$fcs.files[[input$t_6_set_list]])>0)
            {
                update.log(paste("CREATING MUTANT FOR COHORT:", input$t_6_set_list))
                for(i in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[i]])>0 && app.variables$fcs.files[[i]][["set"]]==input$t_6_set_list &&
                       (app.variables$fcs.files[[i]][["type"]]=="CTRL" || app.variables$fcs.files[[i]][["type"]] == "SRC"))
                    {
                        add.mutant(app.variables$fcs.files[[i]][["name"]])
                    }
                }
                update.log("========COHORT CREATED")
            }
        }
            
            
    }, ignoreNULL = T)
    
    observe( #SAVE REDUCTION %
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(length(input$t_6_files_list)>0 && input$t_6_files_list != "")
            {
                if(length(input[["t_6_mutants_list"]])>0 && input[["t_6_mutants_list"]]!="" &&
                   length(app.variables$fcs.files[[input[["t_6_mutants_list"]]]])>0 &&
                   app.variables$fcs.files[[input[["t_6_mutants_list"]]]][["set"]] ==  input[["t_6_set_list"]])
                {
                    current.mutant <- app.variables$loaded.mutant
                    if(!is.null(current.mutant) && length(app.variables$fcs.files[[current.mutant]])>0)
                    {
                        pop.list <- app.variables$populations.list[[input$t_6_files_list]]
                        if(length(pop.list)>0)
                        {
                            app.variables$reduction.percentages[[current.mutant]] <<- sapply(1:length(pop.list), function(j)
                            {
                                return(as.numeric(input[[paste0("t_6_val_",j)]]))
                            })
                        }
                    }
                }
            }
        }
    })
    
    observeEvent(input$t_6_mut_generate, #GENERATE MUTANTS
    {
         shinyjs::disable("t_6_mut_generate")
         if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
         {
             if(!is.na(input$t_6_files_list) && input$t_6_files_list != "")
             {
                 if(length(app.variables$fcs.files[[input$t_6_files_list]]))
                 {
                     mutant.list <- list()
                     for(i in 1:length(app.variables$fcs.files))
                     {
                         if(length(app.variables$fcs.files[[i]])>0)
                         {
                             if(app.variables$fcs.files[[i]][["set"]] ==  input[["t_6_set_list"]] && app.variables$fcs.files[[i]][["type"]] == "MUT")
                             {
                                 mutant.list[[length(mutant.list)+1]] <- names(app.variables$fcs.files)[[i]]
                             }
                         }
                     }
                     if(length(mutant.list)>0)
                     {
                         progress <- Progress$new()
                         progress$set(message = "Generating Mutants", value = 0)
                         lapply(1:length(mutant.list), function(i)
                         {
                             if(length(app.variables$reduction.percentages[[mutant.list[[i]]]]) > 0)
                             {
                                 pop.list <- app.variables$populations.list[[input$t_6_files_list]]
                                 if(length(pop.list)>0)
                                 {
                                     nmb.events.moved <- 0
                                     nmb.pop.increased <- 0
                                     #REDUCTION
                                     lapply(1:length(pop.list), function(j)
                                     {
                                         red.coef <- app.variables$reduction.percentages[[mutant.list[[i]]]][j]
                                         if(red.coef>0)
                                         {
                                             removed.events <- reduce.population(pop.list[[j]], red.coef)
                                             #==
                                             app.variables$used.events[[mutant.list[[i]]]][removed.events] <<- F
                                             #==
                                             tmp.pop <- app.variables$populations.list[[mutant.list[[i]]]][[j]]
                                             nmb.events.moved <<- nmb.events.moved+length(removed.events)
                                             app.variables$populations.list[[mutant.list[[i]]]][[j]] <<-
                                                 tmp.pop[!(tmp.pop%in%removed.events)]
                                         }
                                         else
                                         {
                                             nmb.pop.increased <<- nmb.pop.increased+1
                                         }
                                     })
                                     #AUGMENTATION DE LA TAILLE DES AUTRES POPS
                                     if(nmb.pop.increased>0 && nmb.events.moved>0)
                                     {
                                         nmb.events.per.pop <- as.integer(nmb.events.moved/nmb.pop.increased)
                                         pop.col <- app.variables$fcs.files[[mutant.list[[i]]]][["populations_column"]]
                                         lapply(1:length(pop.list), function(j)
                                         {
                                             red.coef <- app.variables$reduction.percentages[[mutant.list[[i]]]][j]
                                             if(red.coef<=0 && length(pop.list[[j]])>0)
                                             {
                                                 pop <- as.matrix(app.variables$output.matrices[[mutant.list[[i]]]][ pop.list[[j]], ])
                                                 if(length(pop.list[[j]])==1)
                                                 {
                                                     pop <- t(pop)
                                                 }
                                                 
                                                 inc.coef <- 100*nmb.events.per.pop/nrow(pop)
                                                 if(inc.coef>(10^(-9)))
                                                 {
                                                     new.mat <- create.pop.from.pop(pop, inc.coef, unused.columns = pop.col, limited=T)
                                                     first.id <- nrow(app.variables$output.matrices[[mutant.list[[i]]]])
                                                     
                                                     app.variables$output.matrices[[mutant.list[[i]]]] <<- rbind(app.variables$output.matrices[[mutant.list[[i]]]],
                                                                                                                 new.mat)
                                                     app.variables$populations.list[[mutant.list[[i]]]][[j]] <<-
                                                         c(unlist(app.variables$populations.list[[mutant.list[[i]]]][[j]]), (first.id+1):(first.id+nrow(new.mat)))
                                                     app.variables$used.events[[mutant.list[[i]]]] <<- c(app.variables$used.events[[mutant.list[[i]]]],
                                                                                                         rep(T,nrow(new.mat)))
                                                 }
                                             }
                                         })
                                     }
                                 }
                             }
                             progress$inc(1/length(mutant.list), detail=paste0(mutant.list[[i]], " generated"))
                         })
                         progress$set(message = "Mutants Generated", value = 1)
                         delay(500, progress$close())
                     }
                 }
             }
         }
         delay(500, shinyjs::enable("t_6_mut_generate"))
    })
    
    
    
    
    
    
    
    
    #==============================================MIX FILES=======================================================
    
    
    
    
    
    
    
    
    observe( #UPDATE SETS LIST
    {
        update.sets.list(7)
    })
    
    observe( #CLEAR UI
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_7")
        {
            if(length(input[["t_7_source_list"]])>0 && input[["t_7_source_list"]] != "" &&
               length(input[["t_7_target_list"]])>0 && input[["t_7_target_list"]] != "")
            {
                source.name <- input[["t_7_source_list"]]
                target.name <- input[["t_7_target_list"]]
                if(length(app.variables$fcs.files[[source.name]])>0 && length(app.variables$fcs.files[[target.name]])>0)
                {
                    shinyjs::show("t_7_left")
                    shinyjs::show("t_7_right")
                    shinyjs::enable("t_7_mix")
                }
                else
                {
                    shinyjs::hide("t_7_left")
                    shinyjs::hide("t_7_right")
                    shinyjs::disable("t_7_mix")
                }
            }
            else
            {
                shinyjs::hide("t_7_left")
                shinyjs::hide("t_7_right")
                shinyjs::disable("t_7_mix")
            }
        }
        else
        {
            shinyjs::hide("t_7_left")
            shinyjs::hide("t_7_right")
            shinyjs::disable("t_7_mix")
        }
    })
    
    observe( #UPDATE SOURCE AND TARGET LIST
    {
        files.list <- list()
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            for(i in 1:length(app.variables$fcs.files))
            {
                if(length(app.variables$fcs.files[[i]])>0 && length(input$t_7_set_list)>0 && input$t_7_set_list!= "" && 
                   app.variables$fcs.files[[i]][["set"]] == input$t_7_set_list)
                {
                    files.list[[length(files.list)+1]] <- names(app.variables$fcs.files)[[i]]
                }
            }
        }
        updateSelectInput(session, "t_7_source_list", "Select Source File", choices=files.list)
        updateSelectInput(session, "t_7_target_list", "Select Target File", choices=files.list)
    })
    
    output$t_7_pop_tab <- renderUI( #GENERATE POP LIST
    {
        pop.ui <- list()
        if(length(input$t_7_set_list)>0 && input$t_7_set_list!="")
        {
            if(length(input$t_7_source_list)>0 && !is.na(input$t_7_source_list))
            {
                if(length(input$t_7_target_list)>0 && !is.na(input$t_7_target_list))
                {
                    source.name <- input$t_7_source_list
                    target.name <- input$t_7_target_list
                    
                    if(length(app.variables$fcs.files[[source.name]])>0 && length(app.variables$fcs.files[[target.name]])>0)
                    {
                        source.pop.list <- list()
                        source.pop.list.size <- list()
                        for(i in 1:length(app.variables$populations.list[[source.name]]))
                        {
                            if(length(app.variables$populations.list[[source.name]])>0)
                            {
                                pop.name <- names(app.variables$populations.list[[source.name]])[i]
                                source.pop.list[[i]] <- pop.name
                                source.pop.list.size[[pop.name]] <- unlist(app.variables$populations.list[[source.name]][[i]])
                            }
                        }
                        source.pop.list <- unlist(source.pop.list)
                        
                        target.pop.list <- list()
                        target.pop.list.size <- list()
                        for(i in 1:length(app.variables$populations.list[[target.name]]))
                        {
                            if(length(app.variables$populations.list[[target.name]])>0)
                            {
                                pop.name <- names(app.variables$populations.list[[target.name]])[i]
                                target.pop.list[[i]] <- pop.name
                                target.pop.list.size[[pop.name]] <- unlist(app.variables$populations.list[[target.name]][[i]])
                            }
                        }
                        target.pop.list <- unlist(target.pop.list)
                        
                        pop.ui <- tagList(
                            fluidRow(
                                column(
                                    width=2,
                                    p("Population Name")
                                ),
                                column(
                                    width=2,
                                    p("Events")
                                ),
                                column(
                                    width=2,
                                    p("Frequency")
                                ),
                                column(
                                    width=2,
                                    p("Copy from pop ?")
                                ),
                                column(
                                    width=2,
                                    p("% of events to copy")
                                ),
                                column(
                                    width=2,
                                    p("Target pop to copy to")
                                )
                            )
                        )
                        
                        pop.ui <- list(pop.ui, lapply(1:length(source.pop.list), function(i)
                        {
                            tmp.ui <- tagList(
                                fluidRow(
                                    column(
                                        width=2,
                                        h6(source.pop.list[[i]])
                                    ),
                                    column(
                                        width=2,
                                        h6(length(source.pop.list.size[[i]]))
                                    ),
                                    column(
                                        width=2,
                                        h6(paste0(trunc(10000*length(source.pop.list.size[[i]])/sum(app.variables$used.events[[source.name]]))/100, " %"))
                                    ),
                                    column(
                                        width=2,
                                        checkboxInput(paste0("t_7_pop_",i,"_cb"), NULL, value=T)
                                    ),
                                    column(
                                        width=2,
                                        numericInput(paste0("t_7_pop_",i,"_perc"), NULL, value=0)
                                    ),
                                    column(
                                        width=2,
                                        selectInput(paste0("t_7_pop_",i,"_targ"), NULL, choices=target.pop.list)
                                    )
                                )
                            )
                            return(tmp.ui)
                        }))
                    }
                }
            }
        }
        return(pop.ui)
    })
    
    observe( #UPDATE MARKERS LIST
    {
        markers.list <- list()
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(length(input$t_7_source_list)>0 && !is.na(input$t_7_source_list))
            {
                if(length(input$t_7_target_list)>0 && !is.na(input$t_7_target_list))
                {
                    source.name <- input$t_7_source_list
                    target.name <- input$t_7_target_list
                    
                    if(length(app.variables$fcs.files[[source.name]])>0 && length(app.variables$fcs.files[[target.name]])>0)
                    {
                        markers.list <- unlist(app.variables$fcs.files[[source.name]][["markers"]])
                        markers.list <- markers.list[markers.list%in%app.variables$fcs.files[[target.name]][["markers"]]]
                        
                        pop.col.source <- app.variables$fcs.files[[source.name]][["populations_column"]]
                        pop.col.target <- app.variables$fcs.files[[target.name]][["populations_column"]]
                        markers.list <- markers.list[-c(pop.col.source,pop.col.target)]
                    }
                }
            }
        }
        updateSelectInput(session, "t_7_m1", "1st Marker", choices = markers.list, selected = markers.list)
        updateSelectInput(session, "t_7_m2", "2nd Marker", choices = markers.list, selected = markers.list)
    })
    
    observe( #UPDATE POPULATIONS LIST
    {
        source.pop.list <- list()
        target.pop.list <- list()
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(length(input$t_7_source_list)>0 && !is.na(input$t_7_source_list))
            {
                if(length(input$t_7_target_list)>0 && !is.na(input$t_7_target_list))
                {
                    source.name <- input$t_7_source_list
                    target.name <- input$t_7_target_list
                    
                    if(length(app.variables$fcs.files[[source.name]])>0 && length(app.variables$fcs.files[[target.name]])>0)
                    {
                        source.pop.list <- list()
                        for(i in 1:length(app.variables$populations.list[[source.name]]))
                        {
                            if(length(app.variables$populations.list[[source.name]][[i]])>0)
                            {
                                source.pop.list <- list(unlist(source.pop.list), names(app.variables$populations.list[[source.name]])[i])
                            }
                        }
                        source.pop.list <- unlist(source.pop.list)
                        
                        target.pop.list <- list()
                        for(i in 1:length(app.variables$populations.list[[target.name]]))
                        {
                            if(length(app.variables$populations.list[[target.name]][[i]])>0)
                            {
                                target.pop.list <- list(unlist(target.pop.list), names(app.variables$populations.list[[target.name]])[i])
                            }
                        }
                        target.pop.list <- unlist(target.pop.list)
                    }
                }
            }
        }
        updateSelectInput(session, "t_7_pops_list", "Source Population", choices = source.pop.list, selected = source.pop.list)
        updateSelectInput(session, "t_7_popt_list", "Target Population", choices = target.pop.list, selected = target.pop.list)
    })
    
    observeEvent(input$t_7_mix,
    {
        shinyjs::disable("t_7_mix")
        source.pop.list <- list()
        target.pop.list <- list()
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(length(input$t_7_source_list)>0 && !is.na(input$t_7_source_list))
            {
                if(length(input$t_7_target_list)>0 && !is.na(input$t_7_target_list))
                {
                    source.name <- input$t_7_source_list
                    target.name <- input$t_7_target_list
                    
                    if(length(app.variables$fcs.files[[source.name]])>0 && length(app.variables$fcs.files[[target.name]])>0)
                    {
                        progress <- Progress$new()
                        progress$set("Mixing files populations", value=0)
                        
                        source.pop.list <- names(app.variables$populations.list[[source.name]])
                        lapply(1:length(source.pop.list), function(current.source.id)
                        {
                            if(input$t_7_all_pop_cb || input[[paste0("t_7_pop_",current.source.id,"_cb")]])
                            {
                                target.pop.id <- which(names(app.variables$populations.list[[target.name]])==
                                                           input[[paste0("t_7_pop_",current.source.id,"_targ")]])
                                if(length(target.pop.id)>0)
                                {
                                    target.pop.id <- target.pop.id[[1]]
                                    #==
                                    tmp.val <- input[[paste0("t_7_pop_",current.source.id,"_perc")]]
                                    extr.perc <- as.numeric(!is.na(tmp.val))*tmp.val + as.numeric(is.na(tmp.val))*0
                                    if(extr.perc>0)
                                    {
                                        matching.parameters <- which(app.variables$fcs.files[[source.name]][["markers"]]%in%
                                                                         app.variables$fcs.files[[target.name]][["markers"]])
                                        if(length(matching.parameters)>0)
                                        {
                                            matching.parameters <- as.numeric(unlist(matching.parameters))
                                            if(length(matching.parameters)==length(app.variables$fcs.files[[target.name]][["markers"]]))
                                            {
                                                current.source.pop.events <- app.variables$populations.list[[source.name]][[current.source.id]]
                                                nmb.events <- length(current.source.pop.events)
                                                selected.events <- sample(current.source.pop.events, as.integer(extr.perc*nmb.events/100))
                                                #==EXTRACTION DE L'INDICE DE LA POP TARGET DANS LA MATRICE D'EXPRESSION
                                                target.file.pop.col <- app.variables$fcs.files[[target.name]][["populations_column"]]
                                                target.pop <- app.variables$populations.list[[target.name]][[target.pop.id]]
                                                target.pop.matrix.id <- app.variables$output.matrices[[target.name]][as.numeric(unlist(target.pop)),target.file.pop.col]
                                                target.pop.matrix.id <- unique(target.pop.matrix.id)
                                                target.mat.nrow <- nrow(app.variables$output.matrices[[target.name]])
                                                #==
                                                app.variables$output.matrices[[target.name]] <<-
                                                    rbind(app.variables$output.matrices[[target.name]],
                                                          app.variables$output.matrices[[source.name]][as.numeric(unlist(selected.events)),matching.parameters])
                                                #==
                                                app.variables$populations.list[[target.name]][[target.pop.id]] <<-
                                                    c(unlist(app.variables$populations.list[[target.name]][[target.pop.id]]),
                                                      (target.mat.nrow+1):(target.mat.nrow+length(selected.events)))
                                                target.pop <- app.variables$populations.list[[target.name]][[target.pop.id]]
                                                #==
                                                app.variables$output.matrices[[target.name]][as.numeric(unlist(target.pop)),target.file.pop.col] <<-
                                                    target.pop.matrix.id
                                            }
                                        }
                                    }
                                    progress$inc(1/length(source.pop.list), 
                                                 detail=paste0(source.pop.list[[current.source.id]], " to ", 
                                                               input[[paste0("t_7_pop_",current.source.id,"_cb")]]))
                                }
                            }
                        })
                        
                        progress$set("Files Mixed", value=1)
                        delay(500, progress$close())
                    }
                }
            }
        }
        delay(500, shinyjs::enable("t_7_mix"))
    })
    
    output$t_7_plot_1 <- renderPlot( #RENDER PLOT
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(length(input$t_7_source_list)>0 && !is.na(input$t_7_source_list))
            {
                f.name <- input[["t_7_source_list"]]
                if(length(app.variables$fcs.files[[f.name]])>0)
                {
                    if(length(input[["t_7_m1"]])>0 && input[["t_7_m1"]] != "")
                    {
                        if(length(input[["t_7_m2"]])>0 && input[["t_7_m2"]] != "")
                        {
                            if(length(input[["t_7_pops_list"]])>0 && input[["t_7_pops_list"]] != "")
                            {
                                fcs.mat <- app.variables$output.matrices[[f.name]]
                                f.id <- which(names(app.variables$fcs.files)==f.name)[[1]]
                                
                                m1.id <- which(app.variables$fcs.files[[f.name]][["markers"]]==input[["t_7_m1"]])[[1]] 
                                m2.id <- which(app.variables$fcs.files[[f.name]][["markers"]]==input[["t_7_m2"]])[[1]]
                                
                                pop.num <- rep("Other",nrow(fcs.mat))
                                pop.events <- as.numeric(unlist(app.variables$populations.list[[f.name]][[as.numeric(input[["t_7_pops_list"]])]]))
                                if(length(pop.events)>0)
                                {
                                    pop.num[pop.events] <- "Selected Population"
                                }
                                else
                                {
                                }
                                tmp.dataframe <- data.frame(P1=unlist(fcs.mat[,as.integer(m1.id)]),
                                                            P2=unlist(fcs.mat[,as.integer(m2.id)]),
                                                            pop=unlist(pop.num))
                                tmp.dataframe$pop <- as.factor(tmp.dataframe$pop)
                                
                                ggplot(tmp.dataframe, aes(x=P1, y=P2, color=pop)) +
                                    geom_point(size = 0.1, stroke = 0, shape = 16) +
                                    xlab(input[["t_7_m1"]]) +
                                    ylab(input[["t_7_m2"]]) +
                                    theme(legend.position = "top") +
                                    xlim(as.numeric(input[["t_2_min_value"]]), as.numeric(input[["t_2_max_value"]])) +
                                    ylim(as.numeric(input[["t_2_min_value"]]), as.numeric(input[["t_2_max_value"]]))
                            }
                        }
                    }
                }
            }
        }
    })
    
    output$t_7_plot_2 <- renderPlot( #RENDER PLOT
    {
            if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
            {
                if(length(input$t_7_target_list)>0 && !is.na(input$t_7_target_list))
                {
                    f.name <- isolate(input[["t_7_target_list"]])
                    if(length(app.variables$fcs.files[[f.name]])>0)
                    {
                        if(length(input[["t_7_m1"]])>0 && input[["t_7_m1"]] != "")
                        {
                            if(length(input[["t_7_m2"]])>0 && input[["t_7_m2"]] != "")
                            {
                                if(length(input[["t_7_popt_list"]])>0 && input[["t_7_popt_list"]] != "")
                                {
                                    fcs.mat <- app.variables$output.matrices[[f.name]]
                                    f.id <- which(names(app.variables$fcs.files)==f.name)[[1]]
                                    
                                    m1.id <- which(app.variables$fcs.files[[f.name]][["markers"]]==input[["t_7_m1"]])[[1]]
                                    m2.id <- which(app.variables$fcs.files[[f.name]][["markers"]]==input[["t_7_m2"]])[[1]]
                                    
                                    pop.num <- rep("Other",nrow(fcs.mat))
                                    pop.events <- as.numeric(unlist(app.variables$populations.list[[f.name]][[as.numeric(input[["t_7_popt_list"]])]]))
                                    if(length(pop.events)>0)
                                    {
                                        pop.num[pop.events] <- "Selected Population"
                                    }
                                    else
                                    {
                                    }
                                    tmp.dataframe <- data.frame(P1=unlist(fcs.mat[,as.integer(m1.id)]),
                                                                P2=unlist(fcs.mat[,as.integer(m2.id)]),
                                                                pop=unlist(pop.num))
                                    tmp.dataframe$pop <- as.factor(tmp.dataframe$pop)
                                    
                                    ggplot(tmp.dataframe, aes(x=P1, y=P2, color=pop)) +
                                        geom_point(size = 0.1, stroke = 0, shape = 16) +
                                        xlab(input[["t_7_m1"]]) +
                                        ylab(input[["t_7_m2"]]) +
                                        theme(legend.position = "top") +
                                        xlim(as.numeric(input[["t_2_min_value"]]), as.numeric(input[["t_2_max_value"]])) +
                                        ylim(as.numeric(input[["t_2_min_value"]]), as.numeric(input[["t_2_max_value"]]))
                                }
                            }
                        }
                    }
                }
            }
        })
    
    
    
    
    
    
    
    #===========================================VISUALIZE MUTANTS=======================================================
    
    
    
    
    
    
    
    observe( #CLEAR UI
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_4_hm")
        {
            if( (length(input[["t_4_hm_file_list"]])>0 && input[["t_4_hm_file_list"]] != "") ||
                (length(input[["t_4_hm_2_file_list"]])>0 && input[["t_4_hm_2_file_list"]] != "") )
            {
                f.name <- input[["t_4_hm_file_list"]]
                f.name.2 <- input[["t_4_hm_2_file_list"]]
                if(length(app.variables$fcs.files[[f.name]])>0 || length(app.variables$fcs.files[[f.name.2]])>0)
                {
                    shinyjs::show("t_4_hm_right")
                }
                else
                {
                    shinyjs::hide("t_4_hm_right")
                }
            }
            else
            {
                shinyjs::hide("t_4_hm_right")
            }
        }
        else
        {
            shinyjs::hide("t_4_hm_right")
        }
    })
    
    observe( #UPDATE SETS LIST
    {
        update.sets.list("4_hm")
    })
    
    observe( #UPDATE FILES LIST
    {
        tab.id <- "4_hm"
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs==paste0("t_",tab.id)) 
        {
            files.list <- list()
            selected.file <- list()
            if(length(input[[paste0("t_",tab.id,"_set_list")]])>0)
            {
                for(i in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[i]])>0)
                    {
                        if(app.variables$fcs.files[[i]][["set"]] ==  input[[paste0("t_",tab.id,"_set_list")]])
                        {
                            files.list[[length(files.list)+1]] <- names(app.variables$fcs.files)[[i]]
                        }
                    }
                }
            }
            selected.file <- files.list
            if(!is.na(input[[paste0("t_",tab.id,"_file_list")]]) && input[[paste0("t_",tab.id,"_file_list")]]!="")
            {
                selected.file <- as.list(unlist(input[[paste0("t_",tab.id,"_file_list")]])[input[[paste0("t_",tab.id,"_file_list")]]%in%files.list])
            }
            updateSelectInput(session, paste0("t_",tab.id,"_file_list"), "Select 1st File", choices=files.list, selected=selected.file)
            updateSelectInput(session, paste0("t_",tab.id,"_2_file_list"), "Select 2nd File", choices=files.list, selected=selected.file)
        }
    })
    
    observe( #UPDATE MARKERS LIST
    {
        markers.list <- list()
        selected.markers <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_4_hm")
        {
            if(length(input[["t_4_hm_file_list"]])>0 && input[["t_4_hm_file_list"]] != "")
            {
                f.name <- input[["t_4_hm_file_list"]]
                if(length(app.variables$fcs.files[[f.name]])>0)
                {
                    if(app.variables$fcs.files[[f.name]][["set"]] ==  input[["t_4_hm_set_list"]])
                    {
                        markers.list <- app.variables$fcs.files[[f.name]][["markers"]]
                    }
                }
            }
            if(length(input[["t_4_hm_2_file_list"]])>0 && input[["t_4_hm_2_file_list"]] != "")
            {
                f.name <- input[["t_4_hm_2_file_list"]]
                if(length(app.variables$fcs.files[[f.name]])>0)
                {
                    if(app.variables$fcs.files[[f.name]][["set"]] ==  input[["t_4_hm_set_list"]])
                    {
                        if(length(markers.list)>0)
                        {
                            cross.markers <- which(markers.list%in%app.variables$fcs.files[[f.name]][["markers"]])
                            markers.list <- markers.list[unlist(cross.markers)]
                        }
                        else
                        {
                            markers.list <- app.variables$fcs.files[[f.name]][["markers"]]
                        }
                    }
                }
            }
            selected.markers <- markers.list
            if(!is.na(input$t_4_hm_markers_list) && length(input$t_4_hm_markers_list)>0)
            {
                selected.markers <- as.list(unlist(input$t_4_hm_markers_list)[input$t_4_hm_markers_list%in%markers.list])
            }
        }
        updateSelectInput(session, "t_4_hm_markers_list", "Select Markers", choices = markers.list, selected = selected.markers)
    })
    
    output$t_4_hm_sb <- renderSunburst(
    {
        file.plot <- NULL
        if(length(app.variables$fcs.files)>0 && input$tabs=="t_4_hm")
        {
            if(length(input$t_4_hm_file_list)>0 && input$t_4_hm_file_list != "")
            {
                f.name <- input$t_4_hm_file_list
                if(length(app.variables$fcs.files[[f.name]]) > 0)
                {
                    pop.ids <- app.variables$populations.list[[f.name]]
                    
                    tmp.mat <- matrix(ncol = 2, nrow=length(pop.ids))
                    tmp.mat[,1] <- names(pop.ids)
                    tmp.mat[,2] <- sapply(1:length(pop.ids), function(i){return(length(pop.ids[[i]]))})
                    colnames(tmp.mat) <- c("level1", "size")
                    
                    tmp.df <- data.frame(tmp.mat, stringsAsFactors = T)
                    tmp.tree <- d3_nest(tmp.df, value_cols = "size")
                    file.plot <- sunburst(tmp.tree, legend = TRUE)
                }
            }
        }
        return(file.plot)
    })
    
    output$t_4_hm_sb_2 <- renderSunburst(
    {
        file.plot <- NULL
        if(length(app.variables$fcs.files)>0 && input$tabs=="t_4_hm")
        {
            if(length(input$t_4_hm_2_file_list)>0 && input$t_4_hm_2_file_list != "")
            {
                if(length(app.variables$fcs.files[[input$t_4_hm_2_file_list]]) > 0)
                {
                    pop.ids <- app.variables$populations.list[[input$t_4_hm_2_file_list]]
                    tmp.mat <- matrix(ncol = 2, nrow=length(pop.ids))
                    tmp.mat[,1] <- names(pop.ids)
                    tmp.mat[,2] <- sapply(1:length(pop.ids), function(i){return(length(pop.ids[[i]]))})
                    colnames(tmp.mat) <- c("level1", "size")
                    
                    tmp.df <- data.frame(tmp.mat, stringsAsFactors = T)
                    tmp.tree <- d3_nest(tmp.df, value_cols = "size")
                    file.plot <- sunburst(tmp.tree, legend = TRUE)
                }
            }
        }
        return(file.plot)
    })
    
    output$t_4_hm_hm <- renderPlotly(
    {
        file.plot <- NULL
        if(length(app.variables$fcs.files)>0 && input$tabs=="t_4_hm")
        {
            if(length(input$t_4_hm_file_list)>0 && input$t_4_hm_file_list != "")
            {
                f.name <- input$t_4_hm_file_list
                if(length(app.variables$fcs.files[[f.name]]) > 0)
                {
                    pop.ids <- app.variables$populations.list[[f.name]]
                    exp.mat <- app.variables$output.matrices[[f.name]]
                    col.id <- app.variables$fcs.files[[f.name]][["populations_column"]]
                    viewed.markers <- as.numeric(unlist(which(app.variables$fcs.files[[f.name]][["markers"]]%in%
                                                                  input$t_4_hm_markers_list)))
                    viewed.markers <- viewed.markers[viewed.markers!=col.id]
                    
                    tmp.val <- create.pop.table.from.fcs.matrix(exp.mat, pop.ids, viewed.markers, col.id)
                    if(!is.null(tmp.val[[2]]))
                    {
                        file.plot <- heatmaply(tmp.val[[2]], Rowv = T, Colv="Rowv", dendrogram = "none")
                    }
                    
                }
            }
        }
        return(file.plot)
    })
    
    output$t_4_hm_hm_2 <- renderPlotly(
    {
        file.plot <- NULL
        if(length(app.variables$fcs.files)>0 && input$tabs=="t_4_hm")
        {
            if(length(input$t_4_hm_2_file_list)>0 && input$t_4_hm_2_file_list != "")
            {
                f.name <- input$t_4_hm_2_file_list
                if(length(app.variables$fcs.files[[f.name]]) > 0)
                {
                    pop.ids <- app.variables$populations.list[[input$t_4_hm_2_file_list]]
                    exp.mat <- app.variables$output.matrices[[input$t_4_hm_2_file_list]]
                    col.id <- app.variables$fcs.files[[input$t_4_hm_2_file_list]][["populations_column"]]
                    viewed.markers <- as.numeric(unlist(which(app.variables$fcs.files[[input$t_4_hm_2_file_list]][["markers"]]%in%
                                                                  input$t_4_hm_markers_list)))
                    viewed.markers <- viewed.markers[viewed.markers!=col.id]
                    
                    tmp.val <- create.pop.table.from.fcs.matrix(exp.mat, pop.ids, viewed.markers, col.id)
                    if(!is.null(tmp.val[[2]]))
                    {
                        file.plot <- heatmaply(tmp.val[[2]], Rowv = T, Colv="Rowv", dendrogram = "none")
                    }
                }
            }
        }
        return(file.plot)
    })
    
    
    
    
    
    
    
    #=============================================DENSITY PLOTS=========================================================
    
    
    
    
    
    
    
    observe( #CLEAR UI
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_4_hist")
        {
            if(length(input[["t_4_hist_file_list"]])>0 && input[["t_4_hist_file_list"]] != "")
            {
                f.name <- input[["t_4_hist_file_list"]]
                if(length(app.variables$fcs.files[[f.name]])>0)
                {
                    shinyjs::show("t_4_hist_right")
                }
                else
                {
                    shinyjs::hide("t_4_hist_right")
                }
            }
            else
            {
                shinyjs::hide("t_4_hist_right")
            }
        }
        else
        {
            shinyjs::hide("t_4_hist_right")
        }
    })
    
    observe( #UPDATE SETS LIST
    {
        update.sets.list("4_hist")
    })
    
    observe( #UPDATE FILES LIST
    {
        files.list <- list()
        selected.file <- NULL
        selected.file.2 <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(!is.na(input[["t_4_hist_set_list"]])  && input[["t_4_hist_set_list"]]!="")
            {
                for(i in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[i]])>0 && app.variables$fcs.files[[i]][["set"]] ==  input[["t_4_hist_set_list"]])
                    {
                        files.list[[length(files.list)+1]] <- names(app.variables$fcs.files)[[i]]
                    }
                }
            }
            selected.file <- files.list
            if(!is.na(input$t_4_hist_file_list) && input$t_4_hist_file_list!="")
            {
                selected.file <- input$t_4_hist_file_list
            }
            selected.file.2 <- files.list
            if(!is.na(input$t_4_hist_2_file_list) && input$t_4_hist_2_file_list!="")
            {
                selected.file.2 <- input$t_4_hist_2_file_list
            }
        }
        updateSelectInput(session, "t_4_hist_file_list", "Select 1st File", choices=files.list, selected=selected.file)
        updateSelectInput(session, "t_4_hist_2_file_list", "Select 2nd File", choices=files.list, selected=selected.file.2)
    })
    
    output$t_4_hist_plots <- renderUI(
    { 
        plots.ui <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0)
        {
            if(!is.na(input$t_4_hist_file_list)>0 && input$t_4_hist_file_list != "")
            {
                if(length(app.variables$fcs.files[[input$t_4_hist_file_list]])>0)
                {
                    mat <- app.variables$output.matrices[[input$t_4_hist_file_list]]
                    used.events <- app.variables$used.events[[input$t_4_hist_file_list]]
                    plots.ui <- lapply(1:ncol(mat), function(i)
                    {
                        tmp.ui <- NULL
                        tmp.plot <- NULL
                        x <- data.frame(val=mat[used.events,i], group=1)
                        #==
                        if(length(input$t_4_hist_2_file_list)>0 && input$t_4_hist_2_file_list != "" &&
                           length(app.variables$fcs.files[[input$t_4_hist_file_list]])>0)
                        {
                            mat.mut <- app.variables$output.matrices[[input$t_4_hist_2_file_list]]
                            x.2 <- data.frame(val=mat.mut[used.events,i], group=2)
                            x <- rbind(x,x.2)
                            x$group <- as.factor(x$group)
                        }
                        #==
                        tmp.plot <- ggplot(x,aes(x=val, fill=group, group=group, col=group)) + 
                            geom_density(alpha=0.25) +
                            xlab(app.variables$fcs.files[[input$t_4_hist_file_list]][["markers"]][[i]])
                        
                        output[[paste0("t_7_hist_plot_",i)]] <- renderPlot(tmp.plot)
                        
                        tmp.ui <- tagList(
                            column
                            (
                                width=4,
                                plotOutput(paste0("t_7_hist_plot_",i))
                            )
                        )
                        return(tmp.ui)
                    })
                }
            }
        }
        return(plots.ui)
    })
    
    
    
    
    
    
    
    #================================================JOYPLOTS===========================================================
    
    
    
    
    
    
    
    observe( #CLEAR UI
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_4_jp")
        {
            if(length(input[["t_4_jp_set_list"]])>0 && input[["t_4_jp_set_list"]] != "")
            {
                shinyjs::show("t_4_hist_right")
            }
            else
            {
                shinyjs::hide("t_4_hist_right")
            }
        }
        else
        {
            shinyjs::hide("t_4_hist_right")
        }
    })
    
    observe( #UPDATE SETS LIST
    {
        update.sets.list("4_jp")
    })
    
    observe( #UPDATE MARKERS LIST
    {
        markers.list <- list()
        selected.markers <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_4_jp")
        {
            if(!is.na(input$t_4_jp_set_list)>0 && input$t_4_jp_set_list != "")
            {
                for(i in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[i]])>0)
                    {
                        if(app.variables$fcs.files[[i]][["set"]] ==  input[["t_4_jp_set_list"]])
                        {
                            if(length(markers.list)==0)
                            {
                                markers.list <- app.variables$fcs.files[[i]][["markers"]]
                            }
                            cross.markers <- which(markers.list%in%app.variables$fcs.files[[i]][["markers"]])
                            markers.list <- markers.list[unlist(cross.markers)]
                        }
                    }
                }
            }
            selected.markers <- markers.list
            if(!is.na(input$t_4_jp_markers_list) && length(input$t_4_jp_markers_list)>0 && input$t_4_jp_markers_list!="")
            {
                selected.markers <- as.list(unlist(input$t_4_jp_markers_list)[input$t_4_jp_markers_list%in%markers.list])
            }
        }
        updateSelectInput(session, "t_4_jp_markers_list", "Select a Marker", choices = markers.list, selected = selected.markers)
    })
    
    output$t_4_jp_plots <- renderPlot(
    { 
        plot.ui <- NULL
        x <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_4_jp")
        {
            if(!is.na(input$t_4_jp_set_list)>0 && input$t_4_jp_set_list != "" && 
               !is.na(input$t_4_jp_markers_list)>0 && input$t_4_jp_markers_list != "")
            {
                for(i in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[i]])>0 && app.variables$fcs.files[[i]][["set"]] ==  input[["t_4_jp_set_list"]])
                    {
                        used.marker <- unlist(which(app.variables$fcs.files[[i]][["markers"]]==input$t_4_jp_markers_list))
                        used.events <- app.variables$used.events[[i]]
                        if(is.null(x))
                        {
                            x <- data.frame(par=app.variables$output.matrices[[i]][used.events,used.marker], 
                                            file=app.variables$fcs.files[[i]][["name"]])
                        }
                        else
                        {
                            x <- rbind(x, data.frame(par=app.variables$output.matrices[[i]][used.events,used.marker], 
                                                     file=app.variables$fcs.files[[i]][["name"]]))
                        }
                    }
                }
                if(!is.null(x))
                {
                    x$file <- as.factor(x$file)
                    plot.ui <- ggplot(x, aes(x=par,y=file)) + geom_density_ridges2() + 
                        xlab(input$t_4_jp_markers_list) + ylab("density")
                }
            }
        }
        return(plot.ui)
    })
    
    
    
    
    
    
    
    #=============================================DOWNLAOD FILES========================================================
    
    
    
    
    
    
    
    observe( #CLEAR UI
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_5")
        {
            if(length(input[["t_5_set_list"]])>0 && input[["t_5_set_list"]] != "")
            {
                nmb.files <- 0
                for(i in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[i]])>0 && app.variables$fcs.files[[i]][["set"]] == input[["t_5_set_list"]])
                    {
                        nmb.files <- nmb.files+1
                    }
                }
                if(nmb.files>0)
                {
                    shinyjs::enable("t_5_dl")
                }
                else
                {
                    shinyjs::disable("t_5_dl")
                }
            }
            else
            {
                shinyjs::disable("t_5_dl")
            }
        }
        else
        {
            shinyjs::disable("t_5_dl")
        }
    })
    
    observe( #UPDATE SETS LIST
    {
        update.sets.list("5")
    })
    
    output$t_5_ctrl_list <- renderUI(
    {
        ctrl.ui <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_5") 
        {
            if(length(input$t_5_set_list)>0)
            {
                ctrl.ui <- lapply(1:length(app.variables$fcs.files), function(i)
                {
                    tmp.ui <- NULL
                    if(length(app.variables$fcs.files[[i]])>0)
                    {
                        if(app.variables$fcs.files[[i]][["set"]] ==  input$t_5_set_list &&
                           app.variables$fcs.files[[i]][["type"]] == "CTRL")
                        {
                            tmp.ui <- tagList(
                                fluidRow
                                (
                                    style="margin-top:1.7vh",
                                    column(
                                        width=5,
                                        h6(app.variables$fcs.files[[i]][["name"]])
                                    ),
                                    column(
                                        width=2,
                                        h6(paste0("Events: ", nrow(app.variables$output.matrices[[i]])))
                                    ),
                                    column(
                                        width=2,
                                        h6(paste0("Dimensions: ", ncol(app.variables$output.matrices[[i]])))
                                    ),
                                    column(
                                        width=2,
                                        h6(paste0("Populations: ", length(app.variables$populations.list[[i]])))
                                    ),
                                    column(
                                        width=1,
                                        checkboxInput(paste0("t_5_",i,"_cb"), NULL)
                                    )
                                )
                            )
                        }
                    }
                    return(tmp.ui)
                })
            }
        }
        return(ctrl.ui)
    })
    
    output$t_5_mut_list <- renderUI(
    {
        mut.ui <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_5") 
        {
            if(length(input$t_5_set_list)>0)
            {
                mut.ui <- lapply(1:length(app.variables$fcs.files), function(i)
                {
                    tmp.ui <- NULL
                    if(length(app.variables$fcs.files[[i]])>0)
                    {
                        if(app.variables$fcs.files[[i]][["set"]] ==  input$t_5_set_list)
                        {
                            if(app.variables$fcs.files[[i]][["name"]] == app.variables$fcs.files[[i]][["source_ctrl"]] ||
                               app.variables$fcs.files[[i]][["type"]] == "CTRL")
                            {
                                list.mut <- c()
                                for(j in 1:length(app.variables$fcs.files))
                                {
                                    if(length(app.variables$fcs.files[[j]])>0
                                       && app.variables$fcs.files[[j]][["type"]] == "MUT"
                                       && app.variables$fcs.files[[i]][["name"]] == app.variables$fcs.files[[j]][["source_ctrl"]]
                                       && app.variables$fcs.files[[j]][["set"]] == input$t_5_set_list)
                                    {
                                        list.mut <- c(unlist(list.mut),j)
                                    }
                                }
                                if(length(list.mut)>0)
                                {
                                    tmp.ui <- tagList(
                                        fluidRow
                                        (
                                            style="margin-bottom:1.7vh;background-color:#00a65a;color:white;margin-left:0.2%;margin-right:0.2%",
                                            column(
                                                width=1,
                                                h6("CONTROL:")
                                            ),
                                            column(
                                                width=11,
                                                h6(app.variables$fcs.files[[i]][["name"]])
                                            )
                                        )
                                    )
                                    tmp.mut <- lapply(list.mut, function(j)
                                    {
                                        val <- tagList(
                                            fluidRow
                                            (
                                                style="margin-bottom:1.7vh;margin-left:0.2%;margin-right:0.2%",
                                                column(
                                                    width=5,
                                                    h6(app.variables$fcs.files[[j]][["name"]])
                                                ),
                                                column(
                                                    width=2,
                                                    h6(paste0("Events: ", nrow(app.variables$output.matrices[[j]])))
                                                ),
                                                column(
                                                    width=2,
                                                    h6(paste0("Dimensions: ", ncol(app.variables$output.matrices[[j]])))
                                                ),
                                                column(
                                                    width=2,
                                                    h6(paste0("Populations: ", length(app.variables$populations.list[[j]])))
                                                ),
                                                column(
                                                    width=1,
                                                    checkboxInput(paste0("t_5_",j,"_cb"), NULL)
                                                )
                                            )
                                        )
                                        return(val)
                                    })
                                    tmp.ui <- list(tmp.ui, tmp.mut)
                                }
                            }
                        }
                    }
                    return(tmp.ui)
                })
            }
        }
        return(mut.ui)
    })
    
    observeEvent(input$t_5_all, #SELECT ALL
    {
        if(input$t_5_all)
        {
            updateCheckboxInput(session, "t_5_mut_all","Select all", value=T)
            updateCheckboxInput(session, "t_5_ctrl_all","Select all", value=T)
            shinyjs::disable("t_5_mut_all")
            shinyjs::disable("t_5_ctrl_all")
        }
        else
        {
            shinyjs::enable("t_5_mut_all")
            shinyjs::enable("t_5_ctrl_all")
        }
    })
    
    observeEvent(input$t_5_ctrl_all, #SELECT ALL CTRL
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_5") 
        {
            if(length(input$t_5_set_list)>0)
            {
                lapply(1:length(app.variables$fcs.files), function(i)
                {
                    if(length(app.variables$fcs.files[[i]])>0)
                    {
                        if(app.variables$fcs.files[[i]][["set"]] ==  input$t_5_set_list &&
                           app.variables$fcs.files[[i]][["type"]] == "CTRL")
                        {
                            if(input$t_5_ctrl_all)
                            {
                                updateCheckboxInput(session, paste0("t_5_",i,"_cb"), NULL, value=T)
                                shinyjs::disable(paste0("t_5_",i,"_cb"))
                            }
                            else
                            {
                                shinyjs::enable(paste0("t_5_",i,"_cb"))
                            }
                        }
                    }
                })
            }
        }
    })
    
    observeEvent(input$t_5_mut_all, #SELECT ALL MUT
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$t_5_ctrl_all) 
        {
            if(length(input$t_5_set_list)>0)
            {
                lapply(1:length(app.variables$fcs.files), function(i)
                {
                    if(length(app.variables$fcs.files[[i]])>0)
                    {
                        if(app.variables$fcs.files[[i]][["set"]] ==  input$t_5_set_list)
                        {
                            if(app.variables$fcs.files[[i]][["type"]] == "MUT" && 
                               length(app.variables$fcs.files[[app.variables$fcs.files[[i]][["source_ctrl"]]]])>0)
                            {
                                if(input$t_5_mut_all)
                                {
                                    updateCheckboxInput(session, paste0("t_5_",i,"_cb"), NULL, T)
                                    shinyjs::disable(paste0("t_5_",i,"_cb"))
                                }
                                else
                                {
                                    shinyjs::enable(paste0("t_5_",i,"_cb"))
                                }
                            }
                        }
                    }
                })
            }
        }
    })
    
    output$t_5_dl <- downloadHandler(
        filename = function()
        {
            return("output.zip")
        },
        content= function(file)
        {
            shinyjs::disable("t_5_dl")
            files.names <- c()
            nmb.files <- 0
            
            progress <- Progress$new()
            progress$set("Downloading files", value=0)
            if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0) 
            {
                for(i in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[i]])>0 && length(input[[paste0("t_5_",i,"_cb")]])>0 &&
                       !is.na(input[[paste0("t_5_",i,"_cb")]]))
                    {
                        if(input[[paste0("t_5_",i,"_cb")]])
                        {
                            nmb.files <- nmb.files+1
                        }
                    }
                }
            }
            
            progress$inc(1/(nmb.files+2), detail=paste0(nmb.files, " files detected"))
            
            if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0) 
            {
                for(i in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[i]])>0 && length(input[[paste0("t_5_",i,"_cb")]])>0 &&
                       !is.na(input[[paste0("t_5_",i,"_cb")]]))
                    {
                        if(input[[paste0("t_5_",i,"_cb")]])
                        {
                            tmp.dir <- as.character(trunc(as.numeric(Sys.time())))
                            dir.create(tmp.dir)
                            dir.create(paste0(tmp.dir,"/CTRL"))
                            dir.create(paste0(tmp.dir,"/MUT"))
                            #==
                            fcs <- app.variables$fcs.files[[i]][["file"]]
                            fcs@exprs <- app.variables$output.matrices[[i]][app.variables$used.events[[i]],]
                            tmp.path <- NULL
                            if(app.variables$fcs.files[[i]][["type"]]=="MUT")
                            {
                                tmp.name <- paste0(tmp.dir,"/MUT/",app.variables$fcs.files[[i]][["name"]],".fcs")
                                write.FCS(fcs, tmp.name)
                                files.names <- c(unlist(files.names), tmp.name) 
                            }
                            else
                            {
                                tmp.name <- paste0(tmp.dir,"/CTRL/",app.variables$fcs.files[[i]][["name"]],".fcs")
                                write.FCS(fcs, tmp.name)
                                files.names <- c(unlist(files.names), tmp.name)
                            }
                            progress$inc(1/(nmb.files+2), detail=paste0(app.variables$fcs.files[[i]][["name"]], " added to download list"))
                        }
                    }
                }
            }
            zip(file, files.names)
            file.remove(unlist(files.names))
            progress$inc(1/(nmb.files+2), detail="Zip archive generated")
            progress$set("Zip ready", value=1)
            
            delay(500, progress$close())
            delay(500, shinyjs::enable("t_5_dl"))
        }
    )
    
    
    
    
    
    
    
    
    
    #========================================COMPENSATE AND TRANSFORM===================================================
    
    
    
    
    
    
    
    observe( #CLEAR UI
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_8")
        {
            if(length(input[["t_8_set_list"]])>0 && input[["t_8_set_list"]] != "")
            {
                nmb.files <- 0
                for(i in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[i]])>0 && app.variables$fcs.files[[i]][["set"]] == input[["t_8_set_list"]])
                    {
                        nmb.files <- nmb.files+1
                    }
                }
                if(nmb.files>0)
                {
                    shinyjs::show("t_8_left")
                    shinyjs::enable("t_8_comp")
                    shinyjs::enable("t_8_transf")
                }
                else
                {
                    shinyjs::hide("t_8_left")
                    shinyjs::disable("t_8_comp")
                    shinyjs::disable("t_8_transf")
                }
            }
            else
            {
                shinyjs::hide("t_8_left")
                shinyjs::disable("t_8_comp")
                shinyjs::disable("t_8_transf")
            }
        }
        else
        {
            shinyjs::hide("t_8_left")
            shinyjs::disable("t_8_comp")
            shinyjs::disable("t_8_transf")
        }
    })
    
    observe( #UPDATE SETS LIST
    {
        update.sets.list(8)
    })
    
    observe( #UPDATE MARKERS LIST
    {
        markers.list <- list()
        selected.markers <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_8")
        {
            if(!is.na(input$t_8_set_list)>0 && input$t_8_set_list != "")
            {
                for(i in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[i]])>0)
                    {
                        if(app.variables$fcs.files[[i]][["set"]] ==  input[["t_8_set_list"]])
                        {
                            if(length(markers.list)==0)
                            {
                                markers.list <- app.variables$fcs.files[[i]][["markers"]]
                            }
                            cross.markers <- which(markers.list%in%app.variables$fcs.files[[i]][["markers"]])
                            markers.list <- markers.list[unlist(cross.markers)]
                        }
                    }
                }
            }
            selected.markers <- markers.list
            if(!is.na(input$t_8_markers_list) && length(input$t_8_markers_list)>0 && input$t_8_markers_list!="")
            {
                selected.markers <- as.list(unlist(input$t_8_markers_list)[input$t_8_markers_list%in%markers.list])
            }
        }
        updateSelectInput(session, "t_8_markers_list", "Select a Marker", choices = markers.list, selected = selected.markers)
    })
    
    output$t_8_files <- renderUI(
    {
        files.plot <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_8") 
        {
            if(!is.na(input$t_8_set_list) && input$t_8_set_list!="")
            {
                header.ui <- tagList(
                    column
                    (
                        width=7,
                        p("NAME")
                    ),
                    column
                    (
                        width=2,
                        p("COMPENSATED")
                    ),
                    column
                    (
                        width=2,
                        p("TRANSFORMED")
                    ),
                    column
                    (
                        width=1,
                        p("SELECT?")
                    )
                )
                files.plot <- lapply(1:length(app.variables$fcs.files), function(f.id)
                {
                    tmp.ui <- NULL
                    if(length(app.variables$fcs.files[[f.id]])>0)
                    {
                        if(app.variables$fcs.files[[f.id]][["set"]] == input$t_8_set_list)
                        {
                            f.comp <- app.variables$fcs.files[[f.id]][["comp"]]
                            f.transf <- app.variables$fcs.files[[f.id]][["transf"]]
                            logic.to.text <- c("NO","YES")
                            selected <- T
                            if(!is.null(input[[paste0("t_8_cb_",f.id)]]))
                            {
                                seected <- input[[paste0("t_8_cb_",f.id)]]
                            }
                            tmp.ui <- tagList(
                                column
                                (
                                    width=7,
                                    p(app.variables$fcs.files[[f.id]][["name"]])
                                ),
                                column
                                (
                                    width=2,
                                    p(logic.to.text[f.comp+1])
                                ),
                                column
                                (
                                    width=2,
                                    p(logic.to.text[f.transf+1])
                                ),
                                column
                                (
                                    width=1,
                                    checkboxInput(paste0("t_8_cb_",f.id), NULL, value=selected)
                                )
                            )
                        }
                    }
                    return(tmp.ui)
                })
                files.plot <- Filter(Negate(is.null), files.plot)
                if(!is.null(files.plot))
                {
                    files.plot <- list(header.ui, files.plot)
                }
            }
        }
        return(files.plot)
    })
    
    observeEvent(input$t_8_comp,
    {
        progress <- Progress$new()
        progress$set("COMPENSATING FILES", value=0)
        update.log("COMPENSATING FILES")
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_8") 
        {
            if(!is.na(input$t_8_set_list) && input$t_8_set_list!="")
            {
                nmb.files <- 0
                for(f.id in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[f.id]])>0)
                    {
                        if(app.variables$fcs.files[[f.id]][["set"]] == input$t_8_set_list)
                        {
                            nmb.files <- nmb.files+1
                        }
                    }
                }
                if(nmb.files>0)
                {
                    for(f.id in 1:length(app.variables$fcs.files))
                    {
                        if(length(app.variables$fcs.files[[f.id]])>0)
                        {
                            if(app.variables$fcs.files[[f.id]][["set"]] == input$t_8_set_list)
                            {
                                if(!is.null(input[[paste0("t_8_cb_",f.id)]]) && input[[paste0("t_8_cb_",f.id)]])
                                {
                                    tmp.fcs <- app.variables$fcs.files[[f.id]][["file"]]
                                    tmp.fcs <- m.compensate(tmp.fcs)
                                    app.variables$output.matrices[[f.id]] <<- tmp.fcs@exprs
                                    progress$inc(1/nmb.files, detail=paste0(app.variables$fcs.files[[f.id]][["name"]], " compensated"))
                                    update.log(paste("========", app.variables$fcs.files[[f.id]][["name"]],"compensated"))
                                    app.variables$fcs.files[[f.id]][["comp"]] <<- T
                                }
                            }
                        }
                    }
                }
            }
        }
        progress$set("FILES COMPENSATED", value=1)
        update.log("======== FILES COMPENSATED")
        update.log("")
        delay(500, progress$close())
    })
    
    output$t_8_transf_param <- renderUI(
    {
        params.ui <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_8") 
        {
            selected.transform <- as.numeric(input$t_8_transf_type)
            if(selected.transform == 2)
            {
                params.ui <- numericInput("t_8_arcsinh", "Arcsinh cofactor", value="5")
            }
        }
        return(params.ui)
    })
    
    observeEvent(input$t_8_transf,
    {
        progress <- Progress$new()
        progress$set("TRANSFORMING FILES", value=0)
        update.log("TRANSFORMING FILES")
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_8") 
        {
            selected.transform <- as.numeric(input$t_8_transf_type)
            selected.algo <- NULL
            if(selected.transform == 1)
            {
                selected.algo <- m.transform.logicle
                selected.algo.params <- NULL
            }
            else
            {
                selected.algo <- m.transform.asinh
                selected.algo.params <- as.numeric(input$t_8_arcsinh)
            }
            
            if(!is.na(input$t_8_set_list) && input$t_8_set_list!="" && !is.null(selected.algo))
            {
                markers.list <- input$t_8_markers_list
                nmb.files <- 0
                for(f.id in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[f.id]])>0)
                    {
                        if(app.variables$fcs.files[[f.id]][["set"]] == input$t_8_set_list)
                        {
                            nmb.files <- nmb.files+1
                        }
                    }
                }
                if(nmb.files>0)
                {
                    for(f.id in 1:length(app.variables$fcs.files))
                    {
                        if(length(app.variables$fcs.files[[f.id]])>0)
                        {
                            if(app.variables$fcs.files[[f.id]][["set"]] == input$t_8_set_list)
                            {
                                if(!is.null(input[[paste0("t_8_cb_",f.id)]]) && input[[paste0("t_8_cb_",f.id)]])
                                {
                                    tmp.fcs <- app.variables$fcs.files[[f.id]][["file"]]
                                    tmp.fcs <- selected.algo(tmp.fcs, markers.list, selected.algo.params)
                                    app.variables$output.matrices[[f.id]] <<- tmp.fcs@exprs
                                    progress$inc(1/nmb.files, detail=paste0(app.variables$fcs.files[[f.id]][["name"]], " transformed"))
                                    update.log(paste("========", app.variables$fcs.files[[f.id]][["name"]],"transformed"))
                                    app.variables$fcs.files[[f.id]][["transf"]] <<- T
                                }
                            }
                        }
                    }
                }
            }
        }
        progress$set("FILES TRANSFORMED", value=1)
        update.log("======== FILES TRANSFORMED")
        update.log("")
        delay(500, progress$close())
    })
    
    
    
    
    
    
    
    #====================================DECOMPENSATE AND DETRANSFORM===================================================
    
    
    
    
    
    
    
    observe( #CLEAR UI
    {
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_9")
        {
            if(length(input[["t_9_set_list"]])>0 && input[["t_9_set_list"]] != "")
            {
                nmb.files <- 0
                for(i in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[i]])>0 && app.variables$fcs.files[[i]][["set"]] == input[["t_9_set_list"]])
                    {
                        nmb.files <- nmb.files+1
                    }
                }
                if(nmb.files>0)
                {
                    shinyjs::show("t_9_left")
                    shinyjs::enable("t_9_comp")
                    shinyjs::enable("t_9_transf")
                }
                else
                {
                    shinyjs::hide("t_9_left")
                    shinyjs::disable("t_9_comp")
                    shinyjs::disable("t_9_transf")
                }
            }
            else
            {
                shinyjs::hide("t_9_left")
                shinyjs::disable("t_9_comp")
                shinyjs::disable("t_9_transf")
            }
        }
        else
        {
            shinyjs::hide("t_9_left")
            shinyjs::disable("t_9_comp")
            shinyjs::disable("t_9_transf")
        }
    })
    
    observe( #UPDATE SETS LIST
    {
        update.sets.list(9)
    })
    
    observe( #UPDATE MARKERS LIST
    {
        markers.list <- list()
        selected.markers <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_9")
        {
            if(!is.na(input$t_9_set_list)>0 && input$t_9_set_list != "")
            {
                for(i in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[i]])>0)
                    {
                        if(app.variables$fcs.files[[i]][["set"]] ==  input[["t_9_set_list"]])
                        {
                            if(length(markers.list)==0)
                            {
                                markers.list <- app.variables$fcs.files[[i]][["markers"]]
                            }
                            cross.markers <- which(markers.list%in%app.variables$fcs.files[[i]][["markers"]])
                            markers.list <- markers.list[unlist(cross.markers)]
                        }
                    }
                }
            }
            selected.markers <- markers.list
            if(!is.na(input$t_9_markers_list) && length(input$t_9_markers_list)>0 && input$t_9_markers_list!="")
            {
                selected.markers <- as.list(unlist(input$t_9_markers_list)[input$t_9_markers_list%in%markers.list])
            }
        }
        updateSelectInput(session, "t_9_markers_list", "Select a Marker", choices = markers.list, selected = selected.markers)
    })
    
    output$t_9_files <- renderUI(
    {
        files.plot <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_9") 
        {
            if(!is.na(input$t_9_set_list) && input$t_9_set_list!="")
            {
                header.ui <- tagList(
                    column
                    (
                        width=7,
                        p("NAME")
                    ),
                    column
                    (
                        width=2,
                        p("COMPENSATED")
                    ),
                    column
                    (
                        width=2,
                        p("TRANSFORMED")
                    ),
                    column
                    (
                        width=1,
                        p("SELECT?")
                    )
                )
                files.plot <- lapply(1:length(app.variables$fcs.files), function(f.id)
                {
                    tmp.ui <- NULL
                    if(length(app.variables$fcs.files[[f.id]])>0)
                    {
                        if(app.variables$fcs.files[[f.id]][["set"]] == input$t_9_set_list)
                        {
                            f.comp <- app.variables$fcs.files[[f.id]][["comp"]]
                            f.transf <- app.variables$fcs.files[[f.id]][["transf"]]
                            logic.to.text <- c("NO","YES")
                            selected <- T
                            if(!is.null(input[[paste0("t_9_cb_",f.id)]]))
                            {
                                seected <- input[[paste0("t_9_cb_",f.id)]]
                            }
                            tmp.ui <- tagList(
                                column
                                (
                                    width=7,
                                    p(app.variables$fcs.files[[f.id]][["name"]])
                                ),
                                column
                                (
                                    width=2,
                                    p(logic.to.text[f.comp+1])
                                ),
                                column
                                (
                                    width=2,
                                    p(logic.to.text[f.transf+1])
                                ),
                                column
                                (
                                    width=1,
                                    checkboxInput(paste0("t_9_cb_",f.id), NULL, value=selected)
                                )
                            )
                        }
                    }
                    return(tmp.ui)
                })
                files.plot <- Filter(Negate(is.null), files.plot)
                if(!is.null(files.plot))
                {
                    files.plot <- list(header.ui, files.plot)
                }
            }
        }
        return(files.plot)
    })
    
    observeEvent(input$t_9_decomp,
    {
        progress <- Progress$new()
        progress$set("DECOMPENSATING FILES", value=0)
        update.log("DECOMPENSATING FILES")
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_9") 
        {
            if(!is.na(input$t_9_set_list) && input$t_9_set_list!="")
            {
                nmb.files <- 0
                for(f.id in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[f.id]])>0)
                    {
                        if(app.variables$fcs.files[[f.id]][["set"]] == input$t_9_set_list)
                        {
                            nmb.files <- nmb.files+1
                        }
                    }
                }
                if(nmb.files>0)
                {
                    for(f.id in 1:length(app.variables$fcs.files))
                    {
                        if(length(app.variables$fcs.files[[f.id]])>0)
                        {
                            if(app.variables$fcs.files[[f.id]][["set"]] == input$t_9_set_list)
                            {
                                if(!is.null(input[[paste0("t_9_cb_",f.id)]]) && input[[paste0("t_9_cb_",f.id)]])
                                {
                                    tmp.fcs <- app.variables$fcs.files[[f.id]][["file"]]
                                    tmp.fcs <- m.compensate(tmp.fcs)
                                    app.variables$output.matrices[[f.id]] <<- tmp.fcs@exprs
                                    progress$inc(1/nmb.files, detail=paste0(app.variables$fcs.files[[f.id]][["name"]], " decompensated"))
                                    update.log(paste("========", app.variables$fcs.files[[f.id]][["name"]]," decompensated"))
                                    app.variables$fcs.files[[f.id]][["comp"]] <<- F
                                }
                            }
                        }
                    }
                }
            }
        }
        progress$set("FILES DECOMPENSATED", value=1)
        update.log("======== FILES COMPENSATED")
        update.log("")
        delay(500, progress$close())
    })
    
    output$t_9_transf_param <- renderUI(
    {
        params.ui <- NULL
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_9") 
        {
            selected.transform <- as.numeric(input$t_9_transf_type)
            if(selected.transform == 2)
            {
                params.ui <- numericInput("t_9_arcsinh", "Arcsinh cofactor", value="5")
            }
        }
        return(params.ui)
    })
    
    observeEvent(input$t_9_detransf,
    {
        progress <- Progress$new()
        progress$set("DETRANSFORMING FILES", value=0)
        update.log("DETRANSFORMING FILES")
        if(length(app.variables$fcs.files)>0 && length(app.variables$sets.list)>0 && input$tabs=="t_9") 
        {
            selected.transform <- as.numeric(input$t_9_transf_type)
            selected.algo <- NULL
            if(selected.transform == 1)
            {
                selected.algo <- m.transform.logicle
                selected.algo.params <- NULL
            }
            else
            {
                selected.algo <- m.transform.asinh
                selected.algo.params <- as.numeric(input$t_9_arcsinh)
            }
            
            if(!is.na(input$t_9_set_list) && input$t_9_set_list!="" && !is.null(selected.algo))
            {
                markers.list <- input$t_9_markers_list
                nmb.files <- 0
                for(f.id in 1:length(app.variables$fcs.files))
                {
                    if(length(app.variables$fcs.files[[f.id]])>0)
                    {
                        if(app.variables$fcs.files[[f.id]][["set"]] == input$t_9_set_list)
                        {
                            nmb.files <- nmb.files+1
                        }
                    }
                }
                if(nmb.files>0)
                {
                    for(f.id in 1:length(app.variables$fcs.files))
                    {
                        if(length(app.variables$fcs.files[[f.id]])>0)
                        {
                            if(app.variables$fcs.files[[f.id]][["set"]] == input$t_9_set_list)
                            {
                                if(!is.null(input[[paste0("t_9_cb_",f.id)]]) && input[[paste0("t_9_cb_",f.id)]])
                                {
                                    tmp.fcs <- app.variables$fcs.files[[f.id]][["file"]]
                                    tmp.fcs <- selected.algo(tmp.fcs, markers.list, selected.algo.params)
                                    app.variables$output.matrices[[f.id]] <<- tmp.fcs@exprs
                                    progress$inc(1/nmb.files, detail=paste0(app.variables$fcs.files[[f.id]][["name"]], " detransformed"))
                                    update.log(paste("========", app.variables$fcs.files[[f.id]][["name"]]," detransformed"))
                                    app.variables$fcs.files[[f.id]][["transf"]] <<- F
                                }
                            }
                        }
                    }
                }
            }
        }
        progress$set("FILES DETRANSFORMED", value=1)
        update.log("FILES DETRANSFORMED")
        update.log("")
        delay(500, progress$close())
    })
    
    
    
    
    
    
    
    #=================================================LOG===============================================================
    
    
    
    
    
    
    
    output$log <- renderText(
    {
        return(app.variables$log.text)
    })
    
    
    
}
