notebooks <- dir(pattern="*.Rmd")
notebooks <- notebooks[-grep("solution",notebooks)]
data <- c("ozone.csv","cancer.patients.txt","countData.txt","gene.description.txt","gene.expression.txt")

zip("Basic_R_Course.zip", files=c(data,notebooks))