notebooks <- dir(pattern="*.Rmd")
notebooks <- notebooks[-grep("solution",notebooks)]
data <- c("ozone.csv","cancer.patients.txt","countData.txt","gene.description.txt","gene.expression.txt")

file.remove("Basic_R_Course.zip")
zip("Basic_R_Course.zip", files=c(data,notebooks))

solutions <- dir(pattern="solution")
file.remove("Solutions.zip")
zip("Solutions.zip",solutions)
