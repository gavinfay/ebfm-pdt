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
                        tabPanel("Read Me",
                                 includeMarkdown("text_tab.Rmd")),
                        tabPanel("Outcomes", 
                                 #fluidRow(
                                     plotOutput("tsplot"),
                                     h4(textOutput("textout")),
                                     tableOutput("summary1"),
                                 tableOutput("summary2"))#,
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
   
   output$textout <- function() {
     out_table <- report_table(settings = input) 
     
     digits_use <- as.vector(ifelse(out_table[1,]<1,2,0))
     
     if (input$assessType == "stock complex") {
       xx <- out_table %>% 
         filter(bfloor < 1)
       if (nrow(xx)>0) {
         for (i in xx$complex) {
           print(paste0("Stock Complex ",i," is below the biomass floor and a multiplier of ",xx$bfloor[i]," is applied to F."))
         }
       }
       if (nrow(xx)==0) {
         print("Neither stock complex has species assessed to be below the biomass floor, so F is not reduced.")
       }
     }
     
   }
   
   output$summary1 <- function() {
       out_table <- report_table(settings = input) 
       
       digits_use <- as.vector(ifelse(out_table[1,]<1,2,0))

       out_table2 <- out_table
       colnames(out_table2) <- c("Complex","FMSY","BMSY","B_final","F/FMSY","B/BMSY","Catch at FMSY","Floor F multiplier","F","Catch at F","Ceiling","Advice")
       if (input$assessType == "single species")
        colnames(out_table2) <- c("Species","FMSY","BMSY","B_final","F/FMSY","B/BMSY","Catch at FMSY","F","Catch at F","Ceiling","Advice")
       
       out_table2 %>% 
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

   output$summary2 <- function() {
     choices = c("single species", "stock complex")
     params <- as.list(input)
     params$assessType <- choices[!params$assessType == choices]
     
     out_table <- report_table(settings = params) 
     
     digits_use <- as.vector(ifelse(out_table[1,]<1,2,0))
     
     out_table2 <- out_table
     colnames(out_table2) <- c("Complex","FMSY","BMSY","B_final","F/FMSY","B/BMSY","Catch at FMSY","Floor F multiplier","F","Catch at F","Ceiling","Advice")
     if (input$assessType != "single species")
       colnames(out_table2) <- c("Species","FMSY","BMSY","B_final","F/FMSY","B/BMSY","Catch at FMSY","F","Catch at F","Ceiling","Advice")
     
     out_table2 %>% 
       #head(as_tibble(mtcars)) %>% 
       knitr::kable("html", caption = 'Assessment results & catch advice',
                    digits = digits_use, 
                    format.args = list(big.mark = ",", scientific = FALSE)) %>% 
       kableExtra::kable_styling("striped", full_width = F) %>% 
       kableExtra::column_spec(1:ncol(out_table2), color = gray(0.7)) %>% 
       kableExtra::row_spec(0, color = gray(0.7))
                               #bold = TRUE, border_right = TRUE, color = "black", background = "lightgrey")
     # striped = TRUE,  
     # hover = TRUE,
     # caption = 'Assessment results & catch advice')  
     #out_table <- report_table(settings)
     #knitr::kable(out_table, caption = 'Assessment results & catch advice')
   }
   
}

# Run the application 
shinyApp(ui = ui, server = server)

