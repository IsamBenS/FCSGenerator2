# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' @export
FCSGenerator2.run <- function()
{
    library(shiny)
    library(shinydashboard)
    library(shinyjs)
    library(flowCore)
    library(ggplot2)
    library(d3r)
    library(reshape2)
    library(ggridges)
    library(truncnorm)

    library(sunburstR)
    library(heatmaply)
    library(shinyHeatmaply)
    library(tcltk)

    appDir <- system.file("shinyApp", "app", package = "FCSGenerator2")
    if (appDir == "")
    {
        stop("Could not find app directory. Try re-installing `FCSGenerator2`.", call. = FALSE)
    }

    shiny::runApp(appDir, display.mode = "normal", launch.browser = T)
}
