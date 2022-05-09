# function to load the .txt files and remove unnecessary 2nd row 

loadtxt <- function(filepath) {
        file <- read.delim(filepath, skip=2, header=F, sep = "\t")
        filecols <- read.delim(filepath, header=T, sep = "\t")
        colnames(file)<-colnames(filecols)
        return(file)
}
