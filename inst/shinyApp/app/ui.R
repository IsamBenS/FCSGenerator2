ui <- dashboardPage(
    
    dashboardHeader
    (
        title="FCS Generator 2"
    ),
    
    dashboardSidebar
    (
        sidebarMenu
        (
            id="tabs",
            menuItem("Add Files", tabName="t_1"),
            menuItem("Compensate & Transform", tabName="t_8"),
            menuItem("Modify Cohort", tabName="t_cohort",
                     menuSubItem("Generate Mutants", tabName = "t_6"),
                     menuSubItem("Generate Copies", tabName = "t_3")),
            menuItem("Modify Populations", tabName="t_2"),
            menuItem("Mix Files", tabName = "t_7"),
            menuItem("Visualize Populations", tabName = "t_4",
                     menuSubItem("Heatmaps", tabName = "t_4_hm"),
                     menuSubItem("Histograms", tabName = "t_4_hist"),
                     menuSubItem("Joyplots", tabName = "t_4_jp")),
            menuItem("Decompensate & Detransform", tabName="t_9"),
            menuItem("Download", tabName = "t_5"),
            menuItem("Log", tabName = "t_10")
        )
    ),
    
    dashboardBody
    (
        useShinyjs(),
        tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
        fluidRow(
            tabItems
            (
                tabItem
                (
                    tabName="t_1",
                    shinydashboard::box
                    (
                        width=12,title="ADD FILES",solidHeader=T,status="info",
                        shinydashboard::tabBox
                        (
                            width=8,id="t_1_tb",
                            tabPanel
                            (
                                title="Generate New File", value="A",
                                style="padding-right:0",
                                fluidRow
                                (
                                    style="margin-right:0",
                                    column
                                    (
                                        width=2,
                                        numericInput("t_1_nmb_files", "Number of Files", value=1),
                                        numericInput("t_1_nmb_events", "Number of Events", value=1000),
                                        numericInput("t_1_nmb_markers", "Number of Markers", value=2),
                                        numericInput("t_1_nmb_populations", "Number of Populations", value=1)
                                    ),
                                    column
                                    (
                                        width=3,
                                        numericInput("t_1_nmb_rare_populations", "Number of Rare Populations", value=0),
                                        numericInput("t_1_freq_rare_populations", "Total Frequency (Rare Populations)", value=0),
                                        checkboxInput("t_1_transformed", "Transformed (logicle)", value=T),
                                        actionButton("t_1_create", "Generate", width="100%")
                                    ),
                                    column
                                    (
                                        width=7,style="max-height:40vh;overflow:auto;padding:0",
                                        uiOutput("t_1_pop_list")
                                    )
                                )
                            ),
                            tabPanel
                            (
                                title="Import Files",value="B",
                                fluidRow
                                (
                                    column
                                    (
                                        width=4,
                                        actionButton("t_1_select", "Add Files", width="100%")
                                        #actionButton("t_1_compensate", "Compensate Selection", width="100%")
                                    )
                                )
                            )
                        ),
                        # column
                        # (
                        #     width=4, style="max-height:23vh",
                        #     fluidRow
                        #     (
                        #         tableOutput("t_1_fileInfo"),
                        #         style="overflow:auto;max-width:50vw"
                        #     ),
                        #     
                        #     style="overflow:auto;max-width:50vw"
                        # ),
                        column
                        (
                            width=12,
                            uiOutput("t_1_files_main")
                        )
                    )
                ),
                
                tabItem
                (
                    tabName="t_2",style="overflow:auto",
                    column
                    (
                        width=6,
                        shinydashboard::box
                        (
                            width=12,title="OPTIONS",solidHeader=T,status="info",
                            column
                            (
                                width=6, style="max-height:36vh;overflow:auto;",
                                div(
                                    selectInput("t_2_set_list", "Select a Cohort", choices=NULL)
                                ),
                                div(
                                    selectInput("t_2_file_list", "Select a File", choices=NULL)
                                ),
                                div(
                                    selectInput("t_2_markers_list", "Select Markers", choices=NULL, multiple = T),
                                    style="max-height:18vh;overflow:auto"
                                )
                            ),
                            column
                            (
                                width=6,style="max-height:36vh;padding-left:3%",
                                fluidRow
                                (
                                    width=6,
                                    selectInput("t_2_pop_list", "Select a Population", choices=NULL)
                                ),
                                fluidRow
                                (
                                    width=6,
                                    numericInput("t_2_pop_events", "Number of Events", value=100),
                                    column
                                    (
                                        width=4,style="padding:0",
                                        actionButton("t_2_pop_inc", "Add to Pop", width="99%")
                                    ),
                                    column
                                    (
                                        width=4,style="padding:0",
                                        actionButton("t_2_pop_red", "Remove from Pop", width="99%")
                                    ),
                                    column
                                    (
                                        width=4,style="padding:0",
                                        actionButton("t_2_pop_new", "Create New Pop", width="99%")
                                    )
                                )
                            )
                        ),
                        column
                        (
                            width=12,id="t_2_plot_div",
                            shinydashboard::box
                            (
                                width=12,
                                plotOutput("t_2_plot", click = "t_2_plot_click")
                            )
                        )
                    ),
                    column
                    (
                        width=6,id="t_2_right_col",
                        shinydashboard::box
                        (
                            width=12,style="max-height:36vh;",
                            column
                            (
                                width=6,
                                style="overflow:auto;height:33vh;padding-left:10%",
                                h4("Means"),
                                uiOutput("t_2_means")
                            ),
                            column
                            (
                                width=6,
                                style="overflow:auto;height:33vh;padding-left:10%",
                                h4("Standard Deviations"),
                                uiOutput("t_2_sd")
                            )
                        ),
                        shinydashboard::box
                        (
                            width=12,
                            column
                            (
                                width=6,
                                selectInput("t_2_m1", "1st Marker", choices=NULL),
                                selectInput("t_2_m2", "2nd Marker", choices=NULL)
                            ),
                            column
                            (
                                width=6,
                                textInput("t_2_min_value", "Min Value for each Marker", value = -0.5),
                                textInput("t_2_max_value", "Max Value for each Marker", value = 4.5)
                            ),
                            column
                            (
                                width=12,style="overflow:auto",
                                actionButton("t_2_pop_move", "Move Population", width="47%", style="margin-left:3%"),
                                actionButton("t_2_pop_remove", "Delete Population", width="47%")
                            )
                        )
                    )
                ),
                
                tabItem
                (
                    tabName="t_3",
                    column(width=2),
                    shinydashboard::box
                    (
                        width=8,title="GENERATE COPIES",solidHeader=T,status="info",
                        column
                        (
                            width=12, style="height:30vh;padding-top:2%",
                            selectInput("t_3_set_list", "Select a Cohort", choices=NULL),
                            selectInput("t_3_file_list", "Select a File", choices=NULL)
                        ),
                        column
                        (
                            width=12, style="height:30vh",id="t_3_2",
                            column
                            (
                                selectInput("t_3_var_markers_list", "Select Variables Markers", choices=NULL, multiple = T),
                                width = 12, style="display:inline-block;overflow:auto;max-height:11vh"
                            ),
                            column
                            (
                                selectInput("t_3_var_pop_list", "Select Variables Populations", choices=NULL, multiple = T),
                                width = 12, style="display:inline-block;overflow:auto;max-height:9vh"
                            ),
                            div
                            (
                                numericInput("t_3_nmb_ctrl", "Number of Files", value=1),
                                width = "68%", style="display:inline-block;margin-left:1.5%"
                            ),
                            div
                            (
                                actionButton("t_3_control_generate", "Generate", width="100%"),
                                width = "28%", style="display:inline-block;margin-left:1%"
                            )
                        )
                        # column
                        # (
                        #     width=5, style="height:30vh",id="t_3_3"
                        # )
                    ),
                    column(width=2)
                ),
                
                tabItem
                (
                    tabName="t_6",
                    shinydashboard::box
                    (
                        title="OPTIONS",solidHeader=T,status="info",
                        width=3, style="max-height:40vh;padding-top:7%;overflow:auto",
                        selectInput("t_6_set_list", "Select a Cohort", choices=NULL),
                        selectInput("t_6_files_list", "Select a Control File", choices=NULL),
                        actionButton("t_6_mut_add", "Add All Mutants", width="50%", style="margin-left:25%"),
                        actionButton("t_6_mut_generate", "Modify Mutants", width="50%", style="margin-left:25%")
                    ),
                    column
                    (
                        width=5, style="height:auto",
                        uiOutput("t_6_1")
                    )
                ),
                
                tabItem
                (
                    tabName="t_7",
                    column
                    (
                        width=3,
                        shinydashboard::box
                        (
                            width=12,solidHeader=T,status="info",
                            title="FILES SELECTION",
                            selectInput("t_7_set_list", "Select a Cohort", choices=NULL),
                            selectInput("t_7_source_list", "MOVE FROM...", choices=NULL),
                            selectInput("t_7_target_list", "...TO", choices=NULL)
                        ),
                        shinydashboard::box
                        (
                            width=12,solidHeader=T,status="info",
                            title="COPY OPTIONS",
                            checkboxInput("t_7_all_pop_cb", "Copy All Populations", value=T, width = "90%"),
                            numericInput("t_7_all_pop_perc", "Combination %", value = 10),
                            actionButton("t_7_mix", "Mix Files", width="50%", style="margin-left:25%")
                        ),
                        shinydashboard::box
                        (
                            width=12,solidHeader=T,status="info",
                            title="VISUALIZATION OPTIONS",
                            selectInput("t_7_m1", "1st Marker", choices=NULL),
                            selectInput("t_7_m2", "2nd Marker", choices=NULL),
                            selectInput("t_7_pops_list", "Source Population", choices = NULL, selected = NULL),
                            selectInput("t_7_popt_list", "Target Populations", choices = NULL, selected = NULL)
                        )
                    ),
                    shinydashboard::box
                    (
                        width=9,solidHeader=T,status="info",
                        title="Move Populations",
                        column
                        (
                            width=12, style="height:45%;max-height:45vh;overflow:auto",
                            uiOutput("t_7_pop_tab")
                        ),
                        column
                        (
                            width=6,
                            plotOutput("t_7_plot_1")
                        ),
                        column
                        (
                            width=6,
                            plotOutput("t_7_plot_2")
                        )
                    )
                ),
                
                tabItem
                (
                    tabName="t_4_hm",
                    shinydashboard::box
                    (
                        width=3,
                        title="OPTIONS",
                        status="info",
                        solidHeader=T,
                        selectInput("t_4_hm_set_list", "Select a Cohort", choices=NULL),
                        selectInput("t_4_hm_file_list", "Select 1st File", choices=NULL),
                        selectInput("t_4_hm_2_file_list", "Select 2nd File", choices=NULL),
                        selectInput("t_4_hm_markers_list", "Select Markers", choices=NULL, multiple = T)
                    ),
                    column
                    (
                        width=9,id="t_4_hm_right",
                        shinydashboard::box
                        (
                            title="1st File",
                            status="success",
                            solidHeader=T,
                            column
                            (
                                width=12,
                                sunburstOutput("t_4_hm_sb",width = "100%")
                            ),
                            column
                            (
                                width=12,style="margin-top:3vh",
                                plotlyOutput("t_4_hm_hm",width = "100%")
                            )
                        ),
                        shinydashboard::box
                        (
                            title="2nd File",
                            status="danger",
                            solidHeader=T,
                            column
                            (
                                width=12,
                                sunburstOutput("t_4_hm_sb_2",width = "100%")
                            ),
                            column
                            (
                                width=12,style="margin-top:3vh",
                                plotlyOutput("t_4_hm_hm_2",width = "100%")
                            )
                        )
                    )
                ),
                
                tabItem
                (
                    tabName="t_4_hist",
                    shinydashboard::box
                    (
                        width=2,
                        title="OPTIONS",
                        status="info",
                        solidHeader=T,
                        selectInput("t_4_hist_set_list", "Select a Cohort", choices=NULL),
                        selectInput("t_4_hist_file_list", "Select 1st File", choices=NULL),
                        selectInput("t_4_hist_2_file_list", "Select 2nd File", choices=NULL)
                    ),
                    column
                    (
                        width=10,id="t_4_hist_right",
                        shinydashboard::box
                        (
                            title="HISTOGRAMS",
                            status="success",
                            solidHeader=T,
                            width=12,
                            column
                            (
                                width=12,
                                uiOutput("t_4_hist_plots")
                            )
                        )
                    )
                ),
                
                tabItem
                (
                    tabName="t_4_jp",
                    shinydashboard::box
                    (
                        width=2,
                        title="OPTIONS",
                        status="info",
                        solidHeader=T,
                        selectInput("t_4_jp_set_list", "Select a Cohort", choices=NULL),
                        selectInput("t_4_jp_markers_list", "Select a Marker", choices=NULL)
                    ),
                    column
                    (
                        width=10,id="t_4_jp_right",
                        shinydashboard::box
                        (
                            title="Joyplots",
                            status="success",
                            solidHeader=T,
                            width=12,
                            column
                            (
                                width=12,style="min-height:80vh",
                                plotOutput("t_4_jp_plots", height = "80vh")
                            )
                        )
                    )
                ),
                
                tabItem
                (
                    tabName="t_5",
                    column
                    (
                        width=12,
                        shinydashboard::box
                        (
                            width=5,
                            title=div
                            (
                                style="margin:0;padding:0",
                                column
                                (
                                    width=7,
                                    h4("Control Files")
                                ),
                                column
                                (
                                    width=5,
                                    checkboxInput("t_5_ctrl_all","Select all")
                                )
                            ),
                            status="success",
                            solidHeader=T,
                            uiOutput("t_5_ctrl_list")
                        ),
                        shinydashboard::box
                        (
                            width=5,
                            title=div
                            (
                                style="margin:0;padding:0",
                                column
                                (
                                    width=7,
                                    h4("Mutant Files")
                                ),
                                column
                                (
                                    width=5,
                                    checkboxInput("t_5_mut_all","Select all")
                                )
                            ),
                            status="danger",
                            solidHeader=T,
                            uiOutput("t_5_mut_list")
                        ),
                        column
                        (
                            width=2,
                            shinydashboard::box
                            (
                                width=12,solidHeader=T,status="info",
                                title="DOWNLOAD OPTIONS",
                                selectInput("t_5_set_list", "Select a Cohort", choices=NULL),
                                checkboxInput("t_5_all", "Select all"),
                                downloadButton("t_5_dl", "Download Selection")
                            )
                        )
                    )
                ),
                
                tabItem
                (
                    tabName="t_8",
                    column
                    (
                        width=10,id="t_8_left",
                        shinydashboard::box
                        (
                            width=12,
                            uiOutput("t_8_files")
                        )
                    ),
                    shinydashboard::box
                    (
                        width=2,title="Select Files",solidHeader=T,status="info",
                        selectInput("t_8_set_list", "Select a Cohort", choices=NULL),
                        fluidRow
                        (
                            column
                            (
                                width=12,
                                actionButton("t_8_comp", "Compensate", width="100%")
                            )
                        ),
                        fluidRow
                        (
                            style="margin-top:3vh",
                            column
                            (
                                width=12,
                                selectInput("t_8_markers_list", "Markers", choices=NULL, multiple=T)
                            ),
                            column
                            (
                                width=12,
                                selectInput("t_8_transf_type", "Transformation", choices=list("logicle"=1,"arcsinh"=2)),
                                uiOutput("t_8_transf_param"),
                                actionButton("t_8_transf", "Transform", width="100%")
                            )
                        )
                    )
                ),
                
                tabItem
                (
                    tabName="t_9",
                    column
                    (
                        width=10,id="t_9_left",
                        shinydashboard::box
                        (
                            width=12,
                            uiOutput("t_9_files")
                        )
                    ),
                    shinydashboard::box
                    (
                        width=2,title="Select Files",solidHeader=T,status="info",
                        selectInput("t_9_set_list", "Select a Cohort", choices=NULL),
                        fluidRow
                        (
                            column
                            (
                                width=12,
                                actionButton("t_9_decomp", "Decompensate", width="100%")
                            )
                        ),
                        fluidRow
                        (
                            style="margin-top:3vh",
                            column
                            (
                                width=12,
                                selectInput("t_9_markers_list", "Markers", choices=NULL, multiple=T)
                            ),
                            column
                            (
                                width=12,
                                selectInput("t_9_transf_type", "Transformation", choices=list("logicle"=1,"arcsinh"=2)),
                                uiOutput("t_9_transf_param"),
                                actionButton("t_9_detransf", "Detransform", width="100%")
                            )
                        )
                    )
                ),
                
                tabItem
                (
                    tabName="t_10",
                    shinydashboard::box
                    (
                        width=12,title="LOG CONSOLE",solidHeader=T,status="info",
                        verbatimTextOutput("log", placeholder = TRUE)
                    )
                )
            )
        )
        
    )
    
)