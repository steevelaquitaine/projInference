# library(shiny)

# # Define UI for application that draws a histogram
# shinyUI(fluidPage(

#   # Application title
#   titlePanel("Analyse Business"),

#   # # Sidebar with a slider input for the number of bins
#   # sidebarLayout(
#   #   sidebarPanel(
#   #     sliderInput("bins",
#   #                 "Number of bins:",
#   #                 min = 1,
#   #                 max = 50,
#   #                 value = 30)
#   #   ),

#   #   # Show a plot of the generated distribution
#   #   mainPanel(
#   #     plotOutput("distPlot")
#   #   )

#     # Show a summary of the dataset and an HTML table with the
#     # requested number of observations. Note the use of the h4
#     # function to provide an additional header above each output
#     # section.    

 
#   # Sidebar with controls to select a dataset and specify the
#   # number of observations to view. The helpText function is
#   # also used to include clarifying text. Most notably, the
#   # inclusion of a submitButton defers the rendering of output
#   # until the user explicitly clicks the button (rather than
#   # doing it immediately when inputs change). This is useful if
#   # the computations required to render output are inordinately
#   # time-consuming.
#   sidebarLayout(    
    
#     sidebarPanel(
#     ),
          
#     mainPanel(      
#       h4("Observations"),
#       tableOutput("view")  
#     )
#   )
# ))

# Load the ggplot2 package which provides
# the 'mpg' dataset.

#### !!! NEED TO FIX dataset is not transferred between the two files
### it works only when preloaded in working directory !!!!


library('shiny')


# shinyUI(fluidPage(

#   #title
#   titlePanel("censusVis"),
  
#   #sidebar
#   sidebarLayout(
    


#     #### SIDE BAR PANEL #####
#     sidebarPanel(

#       #brief description
#       helpText("Visualize Business Income statement 
#       #   data."),
      
#       #select input to visualize
#       selectInput("var", 
#         label = "Choose what to display",
#         # choices = c("Percent White", "Percent Black",
#         #   "Percent Hispanic", "Percent Asian"),
#         choices = c("all"),
#         # selected = "Percent White"
#         ),
      
#       #slider input
#       sliderInput("range", 
#         label = "Range of interest:",
#         min = 0, max = 100, value = c(0, 100))
#     ),
    



#     #### MAIN PANEL #####
#     mainPanel(

#       #text
#       textOutput("Raw data"),

#       # #Create a new row for the table
#       # fluidRow(
#       #   DT::dataTableOutput("view")
#       # )

#     )
#   ),

# ))


shinyUI(fluidPage(
  titlePanel("censusVis"),
  
  sidebarLayout(
    sidebarPanel(

      #brief description
      helpText("Create demographic maps with 
        information from the 2010 US Census."),

      #select input to visualize
      selectInput("var", 
        label = "Choose what to display",
        # choices = c("Percent White", "Percent Black",
        #   "Percent Hispanic", "Percent Asian"),
        choices = c("all"),
        # selected = "Percent White"
        )
    ),
    

    mainPanel(

      #text
      textOutput("Raw data"),
      
      #Create a new row for the table
      # DT::dataTableOutput("view")
      dataTableOutput("view")
      # fluidRow(
      #   DT::dataTableOutput("view")
      # )

    )
  )
))

# if (!require("DT")) install.packages('DT')

# #to automatically adjust app to browser window dimensions
# fluidPage(
  
#   titlePanel("Business analysis"),


#   sidebarLayout(
    
#     position = "left",
    
#     # sidebarPanel( 
#     #   "sidebar panel"
#     #   selectInput("var", 
#     #     label = "Variables",
#     #     choices = c("All"))
#     #   ),
#     mainPanel(textOutput("Visualize"))
#   ),

#   # # Create a new Row in the UI for selectInputs
#   # fluidRow(
#   #   column(4,
#   #       selectInput("man",
#   #                   "Manufacturer:",
#   #                   c("All",
#   #                     unique(as.character(mpg$manufacturer))))
#   #   ),
#   #   column(4,
#   #       selectInput("trans",
#   #                   "Transmission:",
#   #                   c("All",
#   #                     unique(as.character(mpg$trans))))
#   #   ),
#   #   column(4,
#   #       selectInput("cyl",
#   #                   "Cylinders:",
#   #                   c("All",
#   #                     unique(as.character(mpg$cyl))))
#   #   )
#   # ),

#   # #Create a new Row in the UI for selectInputs
#   fluidRow(
#     # #get unique values for selected dataset's variable
#     column(4,
#         selectInput("sur",
#                     "surname:",
#                     c("All",
#                       unique(as.character(dataset$sur))))
#     )
#   ),
#   #   column(4,
#   #       selectInput("nm",
#   #                   "Name:",
#   #                   c("All",
#   #                     unique(as.character(dataset$Name))))
#   #   ),    
#   #   column(4,
#   #       selectInput("pro",
#   #                   "profession:",
#   #                   c("All",
#   #                     unique(as.character(dataset$profession))))
#   #   )
#   # ),  

#   #Create a new row for the table.
#   fluidRow(
#     DT::dataTableOutput("view")
#   )
# )