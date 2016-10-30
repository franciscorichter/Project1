
library(shiny)
library(ggplot2)
library(P1)
library(ape)
library(apTreeshape)
library(Matrix)
library(parallel)
library(foreach)
library(doParallel)


shinyServer(function(input, output) {

  dataset <- reactive({
    diamonds[sample(nrow(diamonds), input$sampleSize),]
  })
  sim <- eventReactive(input$goButton, {
    phyl2(lambda0=input$lambda,mu0=input$mu,K=(input$lambda-input$mu)/input$beta)$newick
  })

  output$plot <- renderPlot({

    # p <- ggplot(dataset(), aes_string(x=input$x, y=input$y)) + geom_point()
    #
    # if (input$color != 'None')
    #   p <- p + aes_string(color=input$color)
    #
    # facets <- paste(input$facet_row, '~', input$facet_col)
    # if (facets != '. ~ .')
    #   p <- p + facet_grid(facets)
    #

    # if (input$smooth)
    #   p <- p + geom_smooth()
    y = 0
   if(input$var=='lambda'){
     sc = input$scale
      x=seq(0.01,to=sc,by=0.01)
      y = foreach(i = 1:length(x),
              .combine = 'c',
              .multicombine = TRUE) %dopar% {
             llik_st(pars=c(x[i],input$beta,input$mu),setoftrees=S)
      }
      set = data.frame(x,y)
      set = set[!is.na(set$y),]
      p = ggplot(set, aes(x, y))+geom_point()}
    if(input$var=='beta'){
      x=seq(0.001,0.2,by=0.001)
      y = foreach(i = 1:length(x),
              .combine = 'c',
              .multicombine = TRUE) %do% {
            llik_st(pars=c(input$lambda,x[i],input$mu),setoftrees=S)
      }
      set = data.frame(x,y)
      set = set[!is.na(set$y),]
      p = ggplot(set, aes(x, y))+geom_point()}
    if(input$var=='mu'){
      x=seq(0.001,0.5,by=0.01)
      for (i in 1:length(x)){
        y[i] = llik_st(pars=c(input$lambda,input$beta,x[i]),setoftrees=S)
      }
      set = data.frame(x,y)
      set = set[!is.na(set$y),]
      p = ggplot(set, aes(x, y))+geom_point()}
    p <- p +labs(x=input$var,y="log-likelihood") # + ggtitle('llik')
     if (input$smooth)
       p <- p + geom_smooth()
     if (input$jitter)
       p <- p + geom_jitter()

    print(p)

  })

  #s <- sim()
#  ntext <- eventReactive(input$goButton, {
#    input$n
#  })

  output$nText <- renderText({
    ntext()
  })
  #output$plot3 <- renderPlot({
  #  plot(sim())
  #})
  output$plot3 <- renderPlot({
    plot(s$newick)
  })
  output$plot2 <- renderPlot({
    dropex <- drop.fossil(s$newick) # drop extinct species
    plot(dropex)
  })





})
