### R code to generate some of the data and static images used in the materials

## Exercise 4a

png("images/exercise4a.png",width=600,height=300)
par(mfrow=c(1,3))
weather <- read.csv("ozone.csv")
plot(weather$Solar.R, weather$Ozone)
hist(weather$Wind)
boxplot(weather$Ozone ~ weather$Month)
dev.off()

### Image with all the plotting characters

png("images/pch.png")
par(mar=c(0.1,0.1,0.1,0.1))
i <- 0:24

x <- floor(i /5) + 1
y <- i %%5

plot(1:10, type="n", xlim = c(1,5), ylim=c(-1,5),axes=F,xlab="",ylab="")
points(x,y,pch=i+1, cex=2)
text(x,y-0.3,i+1)
dev.off()

## Exercise 4b

png("images/exercise4b.png",width=500,height=500)
par(mfrow=c(2,2))
plot(weather$Solar.R,weather$Ozone, col="orange", pch=16,
     ylab="Ozone level", xlab="Solar Radiation", 
     main="Relationship between ozone level and solar radiation")
hist(weather$Wind, col="purple", xlab="Wind Speed", main="Distribution of Wind Speed", breaks=20, freq=FALSE)
boxplot(weather$Ozone~weather$Month,col=rainbow(5),names=c("May", "Jun", "Jul", "Aug","Sep"),las=2,lab="Ozone Level",main="Distribution of Ozone per-month")
dev.off()

## Exercise 5a

png("images/exercise5a.png",width=900,height=300)
par(mfrow=c(1,3))
plot(weather$Solar.R,weather$Ozone,pch=16,col="lightgreen",ylab="Ozone level",xlab="Solar Radiation")
plot(weather$Wind,weather$Ozone, pch=15,col="steelblue",ylab="Ozone level", xlab="Wind Speed")
plot(weather$Temp,weather$Ozone,pch=17,col="orange", ylab="Ozone level",xlab="Temperature")
dev.off()

## Exercise 5b

png("images/exercise5b.png")
plot(weather$Temp,weather$Ozone, pch=17,
     col="orange", ylab="Ozone level",
     xlab="Temperature")
highO <- which(weather$Ozone > 100)
abline(h=100,lty=2)
points(weather$Temp[highO],weather$Ozone[highO],pch=17,col="red")
dev.off()


## Plot of weights of makes versus females

png("images/males-versus-females.png")
age    <- c(50, 21, 35, 45, 28, 31, 42, 33, 57, 62)
weight <- c(70.8, 67.9, 75.3, 61.9, 72.4, 69.9, 63.5, 
            71.5, 73.2, 64.8)
firstName  <- c("Adam", "Eve", "John", "Mary", "Peter", 
                "Paul", "Joanna", "Matthew", "David", "Sally")
secondName <- c("Jones", "Parker", "Evans", "Davis",
                "Baker","Daniels", "Edwards", "Smith", 
                "Roberts", "Wilson")

consent <- c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE,
             FALSE, TRUE, FALSE, TRUE)

sex <- c("Male", "Female", "Male", "Female", "Male", 
         "Male", "Female", "Male", "Male", "Female")

patients <- data.frame(First_Name = firstName, 
                       Second_Name = secondName, 
                       Full_Name = paste(firstName, secondName), 
                       Sex = factor(sex),
                       Age = age,
                       Weight = weight,
                       Consent = consent,
                       stringsAsFactors = FALSE)

males <- patients$Sex == "Male"
females <- patients$Sex == "Female"
par(mfrow=c(1,2))
plot(patients$Age, patients$Weight,pch=16,type="n")
points(patients$Age[males], patients$Weight[males],pch=16,col="steelblue")
points(patients$Age[females], patients$Weight[females],pch=17,col="orangered1")
legend("topleft", legend=c("M","F"), 
       col=c("steelblue","orangered1"), pch=c(16,17))
boxplot(patients$Weight~patients$Sex)
dev.off()
## Gene Expression dataset


if(!file.exists("gene.expression.txt")){
  
  if(!require(breastCancerNKI) | require(genefilter)) {
    source("http://www.bioconductor.org/biocLite.R")
    biocLite(c("breastCancerNKI","genefilter"))
  }
  data("nki")
  cancer.patients <- pData(nki)[,c("samplename","age","er","grade")]
  genes <- fData(nki)[,c("probe","HUGO.gene.symbol","Cytoband")]
  
  exprs(nki) <- exprs(nki)[!is.na(genes$HUGO.gene.symbol),]
  
  genes <- genes[!is.na(genes$HUGO.gene.symbol),]
  
  ##get the top50 DE genes, plus 500 random
  ps <- NULL
  for(i in 1:nrow(genes)){
    ps[i] <- t.test(exprs(nki)[i,] ~ factor(cancer.patients$er))$p.value
  }
  
  set.seed(070815)
  ind <- order(ps, decreasing = FALSE)[1:50]
  ind <- sort(c(ind, sample(setdiff(1:nrow(genes),ind),500)))
  
  evalues <- exprs(nki)[ind,]
  genes <- genes[ind,]
  library(org.Hs.eg.db) 
  
  chr <- select(org.Hs.eg.db, columns=c("CHR","CHRLOC"),keys = as.character(genes$HUGO.gene.symbol),keytype = "SYMBOL")
  genes$Chromosome <- chr[match(genes$HUGO.gene.symbol, chr[,1]),2]
  genes$Chromosome <- ifelse(!is.na(genes$Chromosome),paste0("chr", genes$Chromosome),NA)
  genes$Start <- abs(chr[match(genes$HUGO.gene.symbol, chr[,1]),3])
  
  genes <- genes[,-3]
  
  final <- !is.na(genes$Chromosome)
  genes <- genes[final,]
  evalues <- evalues[final,]    
  
  
  
  
  
  write.table(evalues, file="gene.expression.txt",quote=FALSE,sep="\t")
  write.table(genes, file="gene.description.txt",quote=FALSE,sep="\t")
  write.table(cancer.patients, file="cancer.patients.txt",quote=FALSE,sep="\t")
}


