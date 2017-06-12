notebooks <- dir(pattern="*.Rmd")
notebooks <- notebooks[-grep("solution",notebooks)]
data <- c("ozone.csv","cancer.patients.txt","countData.txt","patient-info.txt","gene.description.txt","gene.expression.txt","create-patients-data.R")

file.remove("Basic_R_Course.zip")
images <- list.files("images")
zip("Basic_R_Course.zip", files=c(data,notebooks,paste("images/",images,sep="/")))

solutions <- dir(pattern="solution")
file.remove("Solutions.zip")
zip("Solutions.zip",solutions)
