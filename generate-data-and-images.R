### R code to generate some of the data and static images used in the materials

## Get the cleaned patient data from the r-intermediate course

download.file("https://raw.githubusercontent.com/bioinformatics-core-shared-training/r-intermediate/master/patient-data-cleaned.txt",destfile = "patient-data-cleaned.txt")
tmp <- read.delim("patient-data-cleaned.txt")
head(tmp)
set.seed("25052017")
library(dplyr)
rand.ind1 <- sample(1:nrow(tmp),10)
rand.ind2 <- sample(1:nrow(tmp),5)

tmp$Height[rand.ind1] <- NA
tmp$Weight[rand.ind2] <- NA

tmp %>% 
  mutate(Age=floor(runif(nrow(tmp), 23,90))) %>% 
  select(-c(Name, BMI,Overweight,Died,Count,Birth,Date.Entered.Study)) %>% 
write.table("patient-info.txt",quote=FALSE,sep="\t",row.names = FALSE)

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


boxplot(patients$Weight~patients$Sex)
dev.off()

png("images/exercise6.png")
plot(weather$Temp, weather$Ozone, pch=16)
mod1 <- lm(weather$Ozone~weather$Temp)
abline(mod1, col="red", lty=2)
c = coef(mod1)
text(60,150, paste("y = ", round(c[2],2), "x",round(c[1],2),sep=""))
dev.off()

png("images/exercise6b.png")
plot(weather$Temp, weather$Ozone, pch=16)
abline(mod1, col="red", lty=2)
cor = cor(weather$Temp,weather$Ozone,use="c")
cor
text(95,150, paste("r^2 = ", round(cor^2,2)))
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


