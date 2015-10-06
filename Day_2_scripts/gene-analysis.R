### Read the expression data, clinical data, and gene descriptions
genes <- read.delim("gene.description.txt")
subjects <- read.delim("cancer.patients.txt")
evals <- read.delim("gene.expression.txt",stringsAsFactors = FALSE)


## Find out which probe is used to measure ESR1

ind <- match("ESR1", genes[,2])
probe <- genes[ind,1]

## Extract the gene expression measurements for this probe

genevals <- evals[match(probe,rownames(evals)),]

## Produce a boxplot, and do the t-test

boxplot(as.numeric(genevals)~factor(subjects$er))
t.test(as.numeric(genevals)~factor(subjects$er))