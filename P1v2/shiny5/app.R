library(shiny)
sidebarLayout(
  sidebarPanel(
    selectizeInput('main', 'Main title', LETTERS),
    sliderInput('size', 'Point size', min = 0.2, max = 5, value = 1)
  ),
  mainPanel(
    renderPlot(plot(cars, main = input$main, cex = input$size, pch = 19),
               width = 600, height = 400)
  )
)