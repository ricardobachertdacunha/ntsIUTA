#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinydashboard
<<<<<<< HEAD
=======
#' @importFrom shinydashboardPlus dashboardPage dashboardHeader
>>>>>>> 60879cfbc072b46cbbe0945871348c5693aa34bc
#' @import shinyFiles
#' @import shinyWidgets
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # List the first level UI elements here 
<<<<<<< HEAD
    dashboardPage(
      dashboardHeader( 
=======
    shinydashboardPlus::dashboardPage(
      shinydashboardPlus::dashboardHeader( 
>>>>>>> 60879cfbc072b46cbbe0945871348c5693aa34bc
        # title = tagList(span(class = "logo-lg", "ntsIUTA"), 
        # icon("think-peaks")),
        # enable_rightsidebar = FALSE,
        #rightSidebarIcon = "gears",
        dropdownMenuOutput("notifications") # type = "notifications" # the badge status can be info, primary, warning, success,danger
      ),
      dashboardSidebar(
        sidebarMenu(id = "tabs", # Setting id makes input$tabs give the tabName of currently-selected tab
          menuItem("Project Setup", tabName = "project_setup", icon = icon("prescription-bottle")), selected = TRUE,
          menuItem("Make Features", icon = icon("schlix"), startExpanded = TRUE,
            #menuSubItem("MS files", tabName = "setup", icon = icon("prescription-bottle")), selected = TRUE, #vial
            #menuSubItem("Raw Data", tabName = "tic", icon = icon("creative-commons-sampling")),                
            menuSubItem("Peak Picking", tabName = "pp", icon = icon("icicles")), #icicles
            #menuSubItem("QC Check", tabName = "qc", icon = icon("file-contract")), #file-medical-alt
            menuSubItem("Grouping", tabName = "grouping", icon = icon("ioxhost")), #gitter
            menuSubItem("IS Check", tabName = "is", icon = icon("sistrix")), #dashboard
            menuSubItem("Filters", tabName = "filters", icon = icon("eraser")), #tasks #eraser
            menuSubItem("Annotation", tabName = "annotation", icon = icon("tags")), #project-diagram
            menuSubItem("Make Features", tabName = "features", icon = icon("gitter")) #braille #stream
          ),
          menuItem("Workflows", icon = icon("project-diagram"), #tabName = "workflows",
            menuSubItem("Suspect Screening", tabName = "suspects", icon = icon("sistrix")),
            menuSubItem("Trend Analysis", tabName = "trends", icon = icon("think-peaks")),
            menuSubItem("Transformation Products", tabName = "tps", icon = icon("code-branch"))
          ),
          menuItem("Cross Analysis", icon = icon("bar-chart-o"),
            menuSubItem("Instrument Control", tabName = "crossInstrumentControl"),
            menuSubItem("Frequency of Features", tabName = "crossFrequencyFeatures")
          ),          
          menuItem("Applications", icon = icon("th"),
            menuSubItem("Quality Control", tabName = "appQualityControl"),                 
            menuSubItem("Removal Efficiency", tabName = "appRemovalEfficiency")
          )
        )
      ),  
    
      dashboardBody(
        tabItems(
          ## Setup -----
          tabItem("project_setup",
            fluidRow(
              box(title = NULL, solidHeader = TRUE, width = 12,        
                span(
                  actionButton("newproject", label = " New Project", icon = icon("folder-plus"), width = 120),        
                  style = "position:absolute;left: 2.5em;top:0.9em"        
                ),
                span(        
                  shinyDirButton("openproject", icon = icon("folder"), label = " Open Project", style = "width: 120px;",
                    title = "Please select the project folder. It should contain rfiles folder and MS files."
                  ),      
                  style = "position:absolute;left: 13.3em;top:0.9em"          
                ),
                uiOutput("UIsaveprojectbutton"),        
                br(),
                br()          
              ),
              tabBox(side = "left", selected = "Setup", width = 12,
                tabPanel("Setup",
                  uiOutput("UIprojectsetup"),
                  uiOutput("UIsamplelist")
                ),
                tabPanel("Extract",
                  
                ),
                tabPanel("TIC",
                  "Note that when side=right, the tab order is reversed."
                ),
                tabPanel("QC",
                  "Note that when side=right, the tab order is reversed."
                )
                
                
              )
              
              
              
            )       
          )
        )
      ),
      skin = "green"
    )
  )
}

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
 
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'ntsIUTA'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

