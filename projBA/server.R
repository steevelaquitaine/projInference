

#Purpose : Build the app

# install DT package if not yet installed
# if (!require("DT")) install.packages('DT')
library('shiny')
library('readxl')

shinyServer(function(input,output){

  #read SEC (stock exchange commission) excels reports
  # dataset = read_excel('Financial_Report.xls',1)    

  dataset = reactive({

    #if yes, read content
    dataset = read_excel('Financial_Report.xls',1)
    return(dataset)

  })


  # #Filter data based on selections
  # output$view <- DT::renderDataTable(DT::datatable({    
  #   # data <- mpg
  #   # if (input$man != "All") {
  #   #   data <- data[data$manufacturer == input$man,]
  #   # }
  #   # if (input$cyl != "All") {
  #   #   data <- data[data$cyl == input$cyl,]
  #   # }
  #   # if (input$trans != "All") {
  #   #   data <- data[data$trans == input$trans,]
  #   # }
  #   data
  # }))

  # output$tbl = DT::renderDataTable(
  #     iris, options = list(lengthChange = FALSE)
  #   )
  #render dataset and store it in output under the tag "view"
  # output$view = DT::renderDataTable(
  #    dataset, options = list(lengthChange = FALSE)
  #  )

  output$view <- renderDataTable({ dataset })

})