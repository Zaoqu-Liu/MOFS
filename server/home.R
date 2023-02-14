# load data first
load('./db/data/fit.rda')

# get user input
data.input <- reactive({
  req(input$inputFileBtn,cancelOutput = F)
  read.csv(input$inputFileBtn$datapath)
})

# when click to start
observeEvent(input$runBtn,{
  d <- data.input()

  # check data first
  if (ncol(d)>2&colnames(d)[1]!='Feature'&colnames(d)[2]!='Value') {
    # message('')
    shinyWidgets::sendSweetAlert(
      session = session,
      title = "Error",
      type = "error",
      text = "The data you entered does not meet the format, please refer to the InputExample data!" ,
      closeOnClickOutside = TRUE,
      showCloseButton = TRUE
    )
    return()
  }
  if (sum(d$Feature%in%feature$ID)<22) {
    # message('The data you entered does not fit the model!')
    shinyWidgets::sendSweetAlert(
      session = session,
      title = "Error",
      type = "error",
      text = "The data you entered does not fit the model!" ,
      closeOnClickOutside = TRUE,
      showCloseButton = TRUE
    )
    return()
  }
  
  # calculatae 
  d <- d[d$Feature%in%feature$ID,]
  d2 <- as.data.frame(matrix(d$Value,nrow = 1))
  colnames(d2) <- d$Feature
  d2 <- d2[,feature$ID]
  mlppre <- predict(mlpcla,d2)
  mlpprelab <- apply(mlppre , 1, which.max)
  
  # output 
  #message('Probability:\nMOFS1 = ',round(mlppre[1],6),'\nMOFS2 = ',round(mlppre[2],6),'\nMOFS3 = ',round(mlppre[3],6))
  #message('Hence, this sample was identified as MOFS',mlpprelab)
  
  output$resultOutput <- renderUI({
    msg1 <- paste0('Probability:\nMOFS1 = ',round(mlppre[1],3),'; MOFS2 = ',round(mlppre[2],3),'; MOFS3 = ',round(mlppre[3],3))
    msg2 <- paste0('MOFS',mlpprelab)
    
    div(
      p(msg1, style="font-size:90%;"),
      span("Hence, this sample was identified as ",
           style = "color:#000;font-size:100%;font-weight:bold;"),
      
      if(mlpprelab == 1){
        span(msg2,style="font-size:120%;font-weight:bold;background-color:#119da4;color:#FFF")
      }else if(mlpprelab == 2){
        span(msg2,style="font-size:120%;font-weight:bold;background-color:#FF6666;color:#FFF")
      }else{
        # mofs3
        span(msg2,style="font-size:120%;font-weight:bold;background-color:#ffc857;color:#FFF")
      }
    )
  })
  
})

output$demoDownloadBtn <- downloadHandler(
  # filename
  filename = function() {
    "InputExample.csv"
  },
  
  # contect
  content = function(file) {
     data = read.csv("./db/data/InputExample.csv")
     write.csv(data,file = file,row.names = F)
  },
  
  # output file type
  contentType = ".csv"
)  