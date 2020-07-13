#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(deSolve)
library(tidyverse)
library(patchwork)
library(ggtext)
source("functions.R")

ui <- fluidPage(
    
    # Application title
    titlePanel("Demonstration of ecosystem-based catch advice"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            radioButtons("showTimeSeries", "Show data as:",
                         choices = c("time series", "kobe plots"),
                         selected = "time series"),
            radioButtons("useCeiling", "Set catch ceiling?",
                         choices = c("Yes", "No"),
                         selected = "Yes"),
            radioButtons("assessType", "Assessment type?",
                         choices = c("single species", "stock complex"),
                         selected = "stock complex"),
            sliderInput("targetF", "Target F/FMSY",
                        min = 0,
                        max = 1.5,
                        value = 0.75,
                        step = 0.05),
            sliderInput("floorB", "Floor threshold (%max)",
                        min = 0,
                        max = 1,
                        value = 0.2,
                        step = 0.05),
            radioButtons("floorOption", "floor modify with:",
                        choices = c("avg status", "min status"),
                        selected = "avg status")
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            
            # Output: Tabset w/ plot, summary, and table ----
            tabsetPanel(type = "tabs",
                        tabPanel("Outcomes", 
                                 #fluidRow(
                                     plotOutput("tsplot"),
                                     tableOutput("summary"))#,
                        # tabPanel("Biomass", plotOutput("bioplot")),
                        # tabPanel("Catch", plotOutput("catplot")),#,
                        # # tabPanel("Summary", verbatimTextOutput("summary")),
                        #tabPanel("Summary Table", tableOutput("table"))
            )
            #plotOutput("EwEplot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
   output$tsplot <- renderPlot({
       # do_ts_plot(settings = list(
       #     showTimeSeries = input$showTimeSeries,
       #     useCeiling = input$useCeiling,
       #     assessType = input$assessType,
       #     targetF = input$targetF,
       #     floorB = input$floorB,
       #     floorOption = input$floorOption))
      do_ts_plot(settings = input) 
       })
   

   output$summary <- function() {
       out_table <- report_table(settings = input) 
       
       digits_use <- as.vector(ifelse(out_table[1,]<1,2,0))
       
       out_table %>% 
       #head(as_tibble(mtcars)) %>% 
           knitr::kable("html", caption = 'Assessment results & catch advice',
                        digits = digits_use, 
                        format.args = list(big.mark = ",", scientific = FALSE)) %>% 
           kableExtra::kable_styling("striped", full_width = F)
       # striped = TRUE,  
       # hover = TRUE,
       # caption = 'Assessment results & catch advice')  
       #out_table <- report_table(settings)
       #knitr::kable(out_table, caption = 'Assessment results & catch advice')
   }
 
}

# Run the application 
shinyApp(ui = ui, server = server)

