library(shiny)
library(P1)#
library(phytools)
# Only run examples in interactive R sessions
if (interactive()) {
  
  ui <- fluidPage(
    sidebarPanel(
      numericInput("lambda", "Lambda:", min = 0, max = 10, value = 0.8),
      numericInput("K", "K:", min = 1, max = 200, value = 40),
      numericInput("mu", "Mu:", min = 0, max = 2, value = 0.1),
      numericInput("tt", "Crown time:", min = 1, max = 30, value = 15),
      br(),
      actionButton("goButton", "Simulate tree"),
      p("Click the button to update the value displayed in the main panel.")
    ),
    sidebarPanel(
      checkboxInput("drop", "Drop extinct:", FALSE),
      checkboxInput("tip", "Show tip label:", FALSE)
    #  numericInput("tip", "Show tip label:", min = 0, max = 2, value = 0.1),
      #br(),
  #    actionButton("goButton", "Simulate tree"),
   #   p("Click the button to update the value displayed in the main panel.")
    ),
  #  sliderInput("obs", "Number of observations", 0, 1000, 500),
#    actionButton("goButton", "Go!"),
    plotOutput("distPlot"),
    plotOutput('ltt')
  )
  
  server <- function(input, output) {
    randomVals <- eventReactive(input$goButton, {
      seed=round(runif(1,1,10000000))
      phyl2(tt=input$tt, lambda0=input$lambda,mu0=input$mu,K=input$K, seed=seed)
      
    })
    output$distPlot <- renderPlot({
      # Take a dependency on input$goButton. This will run once initially,
      # because the value changes from NULL to 0.
      #input$goButton
      #seed=round(runif(1,1,10000000))
      # Use isolate() to avoid dependency on input$obs
      #dist <- isolate(phyl2(tt=input$tt, lambda0=input$lambda,mu0=input$mu,K=input$K, seed=seed))
  #    par(mfrow=c(2,1))
      if (input$drop){
        dropex <- drop.fossil(randomVals()$newick) # drop extinct species
        plot(dropex, show.tip.label = input$tip)
      }
      else{
        plot(randomVals()$newick,show.tip.label = input$tip)
      }
 #     ltt(randomVals()$newick)
    })
    output$ltt <- renderPlot({
      ltt(randomVals()$newick)
    })
  }
  
  shinyApp(ui, server)
  
}