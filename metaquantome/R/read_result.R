
# read in result table
read_result <- function(file){
    df <- read.delim(file, sep="\t")
    return(df)
}

