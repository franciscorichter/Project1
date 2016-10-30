library(shiny)
library(ggplot2)

dataset <- diamonds

shinyUI(fluidPage(
  
  title = "Parameters",
  
  plotOutput('plot'),
  
  hr(),
  
  fluidRow(
    column(3,
           h4("Parameters"),
           sliderInput("lambda", "lambda",
                       min = 0.001, max = 1, value = 0.8, step = 0.01),
           sliderInput("beta", "beta",
                       min = -1, max = 0.3, value = 0.0175, step = 0.0001),
           sliderInput("mu", "mu",
                       min = 0.01, max = 0.9, value = sum(1-E)/(sum(n*t)), step = 0.01),
           br(),
           checkboxInput('jitter', 'Jitter'),
           checkboxInput('smooth', 'Smooth')
    ),
    column(4, offset = 1,
           selectInput('var', 'X-axe', c('lambda','beta','mu')),
           sliderInput('scale', 'scale', value=1,min=0.5,max=15,step=0.5),
           selectInput('color', 'Color', c('None', names(dataset)))
    ),
    # column(4,
    #        fluidPage(
    #          # tags$head(tags$script(src = "message-handler.js")),
    #          actionButton("go", "Simulate Tree")
    #         
    #          )
    #        ),
    sidebarPanel(
      numericInput("lambda", "Lambda:", min = 0, max = 10, value = 0.8),
      numericInput("beta", "Beta:", min = -1, max = 1, value = 0.0175),
      numericInput("mu", "Mu:", min = 0, max = 2, value = 0.1),
      br(),
      actionButton("goButton", "Simulate tree"),
      p("Click the button to update the value displayed in the main panel.")
    )
           #selectInput('facet_row', 'Facet Row',
          #             c(None='.', names(diamonds[sapply(diamonds, is.factor)]))),
          # selectInput('facet_col', 'Facet Column',
          #             c(None='.', names(diamonds[sapply(diamonds, is.factor)])))
           
    #)
  ),
  
  plotOutput('plot3'),
  
  plotOutput("plot2")
))