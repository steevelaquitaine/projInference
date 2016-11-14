

# getOnlineFile.R
#
#author : steeve laquitaine
#input :
#	url <- "https://www.sec.gov/Archives/edgar/data/320193/000119312515356351/Financial_Report.xlsx"
#   savedPath <-'~/Desktop/projBusinAna_dev/data/is2016.xls'

getOnlineFile <- function(url,savedPath){
	download.file(url,savedPath, mode="wb")    
}




