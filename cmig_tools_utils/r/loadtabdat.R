# function to load the .txt files and remove unnecessary 2nd row 

loadtabdat <- function(filepath,filename) {
  
        allfiles<-dir(filepath)
        if (str_detect(allfiles[1],'txt')==TRUE){
            file <- read.delim(paste0(filepath,'/',filename,'.txt'), skip=2, header=F, sep = "\t")
            filecols <- read.delim(paste0(filepath,'/',filename,'.txt'), header=T, sep = "\t")
            colnames(file)<-colnames(filecols)
        } else if (str_detect(allfiles[1],'csv')==TRUE){
            file<-read.csv(paste0(filepath,'/',filename,'.csv'))
        }
        return(file)
}
