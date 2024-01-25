# Preparation

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# install.packages("RPMM")
# install.packages("sqldf")
# install.packages("WGCNA")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("impute")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("preprocessCore")

# Load necessary packages.

library(WGCNA)
library(sqldf)
library(impute)

# Read additional data files for probe annotation.

probeAnnotation21kdatMethUsed <- read.csv("AdditionalFile22probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k <- read.csv("AdditionalFile21datMiniAnnotation27k.csv")
datClock <- read.csv("AdditionalFile23predictor.csv")

# Read in the DNA methylation data (beta values)

# Load, clean, and preview

dat0 <- read.csv.sql("AdditionalFile26MethylationDataExample55.csv");
nSamples <- dim(dat0)[[2]]-1
nProbes <- dim(dat0)[[1]]
dat0[,1] <- gsub(x = dat0[,1], pattern = "\"", replacement = "")
cat("Number of Samples:", nSamples, "\n")
cat("Number of Probes:", nProbes, "\n")

# Validate data structure

if (file.exists("LogFile.txt"))
  file.remove("LogFile.txt")
file.create("LogFile.txt")

source("ErrorMessages.r")

DoNotProceed <- FALSE

if (nSamples == 0) {
  DoNotProceed <- TRUE
  cat(paste(error_messages$no_samples_error), file=log_file, append=TRUE)
}

if (nProbes == 0) {
  DoNotProceed <- TRUE
  cat(paste(error_messages$no_probes_error), file=log_file, append=TRUE)
}

if (nSamples > nProbes) {
  cat(paste(error_messages$samples_exceed_probes_warning), file=log_file, append=TRUE)
}

if (is.numeric(dat0[,1])) {
  DoNotProceed <- TRUE
  cat(paste(error_messages$numeric_first_column_error, dat0[1:3,1]), file=log_file, append=TRUE)
}

if (!is.character(dat0[,1])) {
  cat(paste(error_messages$non_character_first_column_warning, dat0[1:3,1]), file=log_file, append=TRUE)
}

datout <- data.frame(Error = c(error_messages$input_error), Comment=c("", "email Steve Horvath."))

if (DoNotProceed) {
  stop("DoNotProceed")
}

nonNumericColumn <- rep(FALSE, dim(dat0)[[2]]-1)
for (i in 2:dim(dat0)[[2]] ) {
  nonNumericColumn[i-1]=! is.numeric(dat0[,i])
}
if (sum(nonNumericColumn) >0) {
  cat(paste(error_messages$non_numeric_columns, colnames(dat0)[-1][ nonNumericColumn], error_messages$non_numeric_columns_hint),file="LogFile.txt",append=TRUE)
}

XchromosomalCpGs <- as.character(probeAnnotation27k$Name[probeAnnotation27k$Chr=="X"])
selectXchromosome <- is.element(dat0[,1], XchromosomalCpGs)
selectXchromosome[is.na(selectXchromosome)] <- FALSE
if (sum(selectXchromosome) >= 500)  {
  meanXchromosome = rep(NA, dim(dat0)[[2]]-1)
  meanXchromosome = as.numeric(apply( as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE))
}
if (sum(is.na(meanXchromosome)) > 0) {
  cat(paste(error_messages$x_missing_probes),file="LogFile.txt",append=TRUE)  # Output message
}

match1 <- match(probeAnnotation21kdatMethUsed$Name, dat0[,1])
if (sum( is.na(match1))>0) { 
  missingProbes <- probeAnnotation21kdatMethUsed$Name[!is.element( probeAnnotation21kdatMethUsed$Name , dat0[,1])]
  DoNotProceed <- TRUE;
  cat(paste(error_messages$missing_probes_1, length(missingProbes), error_messages$missing_probes_2, paste( missingProbes, sep="",collapse=", ")),file="LogFile.txt",append=TRUE)
}

# Filter only 21k relevant probes and ensure they are numeric

match1 <- match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
if (sum(is.na(match1))>0)
  stop(paste(sum(is.na(match1)), "CpG probes cannot be matched"))

dat1 <- dat0[match1,]

asnumeric1 <- function(x) {
  as.numeric(as.character(x))
}

dat1[,-1] <- apply(as.matrix(dat1[,-1]), 2, asnumeric1)

# Making age prediction on samples

set.seed(1)

source("AdditionalFile24NORMALIZATION.r")

normalizeData <- TRUE

source("AdditionalFile25StepwiseAnalysis.r")

write.table(datout, "Output.csv", row.names = F, sep = ",")

# Explanation of the output

head(datout, n=2)

# Estimate of DNAmAge

signif(datout$DNAmAge,2)

# Normalized DNA methylation data preview

dim(datMethUsedNormalized)
dat0UsedNormalized=data.frame(CpGName=colnames(datMethUsedNormalized), data.frame(t(datMethUsedNormalized) ))
dat0UsedNormalized[1:5,1:5]

# Relate DNAm age to chronological age

datSample=read.csv("AdditionalFile27SampleAnnotationExample55.csv")
head(datSample[, 2:6], n = 3)
DNAmAge <- datout$DNAmAge
medianAbsDev <- function(x, y) median(abs(x - y), na.rm = TRUE)
medianAbsDev1 <- signif(medianAbsDev(DNAmAge, datSample$Age), 2)
par(mfrow = c(1, 1))

verboseScatterplot(
  DNAmAge,
  datSample$Age,
  xlab = "DNAm Age",
  ylab = "Chronological Age",
  main = paste("All, err=", medianAbsDev1)
)

abline(0, 1)

# Comparing predicted gender with known gender

predictedGender
datSample$Female

# Prepare variables for analysis

attach(datSample)
AgeAccelerationDiff = DNAmAge - datSample$Age

restNonMissing = !is.na(DNAmAge) & !is.na(Age)
AgeAccelerationResidual = rep(NA, length(Age))
if (sum(restNonMissing, na.rm = TRUE) > 3) {
  DNAmAge_numeric <- as.numeric(datout$DNAmAge)
  Age_numeric <- as.numeric(datSample$Age)
  fit <- lm(DNAmAge_numeric ~ Age_numeric, subset = restNonMissing)
  residuals_from_fit <- residuals(fit)
  AgeAccelerationResidual[restNonMissing] <- residuals_from_fit
}

DiseaseLabel = ifelse(datSample$diseaseStatus == 1, "A", "C")
DiseaseColor = ifelse(datSample$diseaseStatus == 1, "red", "black")

head(AgeAccelerationDiff)
head(AgeAccelerationResidual)

# Plot A: DNAm Age vs. Chronological Age

verboseScatterplot(DNAmAge, Age, xlab = "DNAm Age", ylab = "Chronological Age", main = paste("A err=", medianAbsDev1), type = "n")
abline(0, 1)
text(DNAmAge, Age, lab = DiseaseLabel, col = DiseaseColor)

# Plot B: Chronological Age vs. Disease Status

verboseBarplot(Age, DiseaseLabel,xlab="Disease Status",main="B Age",ylab="Chronological Age",col=c("red","grey"))

# Plot C: DNAm Age vs. Disease Status

verboseBarplot(DNAmAge, DiseaseLabel,xlab="Disease Status",main="C DNAm Age",ylab="DNAm Age",col=c("red","grey"))

# Plot D: Age Acceleration Difference vs. Disease Status

verboseBarplot(AgeAccelerationDiff, DiseaseLabel,main="D Disease Status", xlab="Disease Status",col=c("red","grey"))

# Plot E: Cause of Death vs. Age Acceleration Difference

verboseBarplot(AgeAccelerationDiff, CauseofDeath, main="E Cause of Death",  xlab="Cause of Death")

# Plot F: Post Mortem Interval vs. Age Acceleration Difference

verboseScatterplot(PostMortemInterval, AgeAccelerationDiff, xlab="Post Mortem Interval",main=paste("F PMI, err=", medianAbsDev1),type="n")
text(PostMortemInterval,  AgeAccelerationDiff, lab= DiseaseLabel, col= DiseaseColor   )

# Plot G: Age Acceleration Residual vs. Disease Status

verboseBarplot(AgeAccelerationResidual, DiseaseLabel,main="G Disease Status", xlab="Disease Status",  col=c("red","grey"))

# Plot H: Age Acceleration Residual vs. Cause of Death

verboseBarplot(AgeAccelerationResidual, CauseofDeath, main="H Cause of Death",  xlab="Cause of Death")

# Plot I: Age Acceleration Residual vs. Post Mortem Interval

verboseScatterplot(PostMortemInterval, AgeAccelerationResidual, xlab="Post Mortem Interval",main=paste("I PMI, err=", medianAbsDev1),type="n")
text(PostMortemInterval, AgeAccelerationResidual, lab=DiseaseLabel, col=DiseaseColor)
