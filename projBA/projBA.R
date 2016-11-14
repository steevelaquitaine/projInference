

# author : steeve laquitaine
#purpose : make web app for business analysis
#
# required packages :
#	- shiny
#   - DT

#set default path
setwd('~/Dropbox/myDropbox/Codes/projBA/')

#---------------------- Dowload file from the internet
url <- "https://www.sec.gov/Archives/edgar/data/320193/000119312515356351/Financial_Report.xlsx"
download.file(url, '~/Dropbox/myDropbox/Codes/projBA/data/is2016.xlsx', mode="wb")

#--------------------- Read multiple excel file with multiple sheets
#read relevant file sheets 
#(cleanup) Omit column names
library(readxl)    
operations <- read_excel("data/is2016.xlsx", sheet = 2, col_names = FALSE,skip = 2)
incomes <- read_excel("data/is2016.xlsx", sheet = 3, col_names = FALSE,skip = 2)
balanceSheet <- read_excel("data/is2016.xlsx", sheet = 5, col_names = FALSE,skip = 2)


#---------------------- Run web app
library(shiny)
library('readxl')
library(DT)
runApp(list(

    #------------
    #building app  
    #------------
    server = function(input,output) {

      #Select key data
      #read excels financial report dataset
      #set header
      #Omit variable column 1
      NetIncome = which(incomes=='Net income')            
      NetSales = which(operations=='Net sales')            
      dataset = incomes[c(NetIncome),-1]            
      dataset[2,] = operations[c(NetSales),-1]            
      rownames(dataset) = c("NetIncome","NetSales")

      #store in output to be sent to ui
      #format dataset  
      output$dataset = renderDataTable({                        
        DT::datatable(dataset, rownames = TRUE,selection = 'none', options = list(dom = 't'))
      })

     #Plot "Net income"
     #transpose table
      output$plotNI <- renderPlot({                
        dataset_t = t(dataset)
        plot(dataset_t[,"NetIncome"],type="l",xlab="date",ylab="Net Income (Million USD)")
        par(new=T)
        plot(dataset_t[,"NetSales"],type="l",col="green")
        par(new=T)        
      })
    },

    #--------------
    #user interface
    #--------------
    #Sidebar
    #side panel
    #main panel
    ui = fluidPage(      
      sidebarLayout(          
          position = c("left","right"),                    
          sidebarPanel(dataTableOutput('dataset')),          
          mainPanel(plotOutput('plotNI'))          
      ) 
    )  
 )
)