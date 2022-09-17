#!/usr/bin/R

##############################################################################
#
# MODULE:       r.avaflow.roc.R
#
# AUTHOR:       Martin Mergili
#
# PURPOSE:      The simulation model for avalanche and debris flows
#               Script for ROC plot
#
# COPYRIGHT:    (c) 2013 - 2021 by the author
#               (c) 2020 - 2021 by the University of Graz
#               (c) 2013 - 2021 by the BOKU University, Vienna
#               (c) 2015 - 2020 by the University of Vienna
#               (c) 1993 - 2021 by the R Development Core Team
#
# VERSION:      20210525 (25 May 2021)
#
#               This program is free software under the GNU General Public
#               License (>=v2).
#
##############################################################################

# Loading libraries
library(rgdal)
library(ROCR)

# Importing arguments (defined in r.avaflow.py)
args <- commandArgs(trailingOnly = TRUE)
temppath <- args[1]
prefix <- args[2]
mode <- as.integer(args[3])
normalization <- as.integer(args[4])
type <- as.integer(args[5])

# Creating plot file
if (normalization == 1 && type == 1) rocplot = paste(prefix, "_results/", prefix, 
    "_plots/", prefix, "_roc_iii.png", sep = "")
if (normalization == 2 && type == 1) rocplot = paste(prefix, "_results/", prefix, 
    "_plots/", prefix, "_roc_iii_n.png", sep = "")
if (normalization == 1 && type == 2) rocplot = paste(prefix, "_results/", prefix, 
    "_plots/", prefix, "_roc_dii.png", sep = "")
if (normalization == 2 && type == 2) rocplot = paste(prefix, "_results/", prefix, 
    "_plots/", prefix, "_roc_dii_n.png", sep = "")
png(filename = rocplot, width = 12, height = 10, units = "cm", res = 300)

# Defining margins
par(mar = c(3.2, 3.2, 0.5, 0.5))

# Initializing ROC plot
plot(x = c(-2, 2), y = c(-2, 2), col = "gray", lty = 3, lwd = 1, xlim = c(-0.05, 
    1), ylim = c(0, 1.02), axes = FALSE, xlab = NA, ylab = NA)
clip(0, 1, 0, 1)
abline(a = 0, b = 1, col = "gray", lty = 2, lwd = 1)
clip(-1, 2, -1, 2)

# Initializing vectors
arocb <- vector("numeric", 1)
aroc <- vector("expression", 1)
lwidth <- vector("character", 1)
lltp <- vector("integer", 1)
lcol <- vector("character", 1)

# Defining formatting of ROC plot
i = 1
mval = 0
ia = ""
lintsp = 1
mstring = ""
lwidth[i] = "2"
lltp[i] = 1
if (mode == 1) lcol[i] = "aquamarine4"
if (mode == 2) lcol[i] = "cornflowerblue"
if (mode == 3) lcol[i] = "darkorange"
if (mode == 4) lcol[i] = "brown1"

# Generating data frames from raster maps
observ <- paste(temppath, "/observed", mstring, ".asc", sep = "")
ind <- paste(temppath, "/index", mstring, ".asc", sep = "")
observed <- readGDAL(fname = observ)  #observation
index <- readGDAL(fname = ind)  #index

# Preprocessing data frames
index$band1[index$band1 < 0] <- NaN
if (normalization == 2) index$band1[index$band1 == 0 & observed$band1 == 0] <- NaN
observed$band1[observed$band1 > 0] <- 1
observed$band1[observed$band1 < 0] <- NaN
index$band1[is.na(observed$band1) == TRUE] <- NaN
observed$band1[is.na(index$band1) == TRUE] <- NaN

indexf <- index$band1[!is.na(index$band1)]
observedf <- observed$band1[!is.na(observed$band1)]

if (normalization == 2) {

    npos <- sum(observedf == 1, na.rm = TRUE)
    nneg <- sum(observedf == 0, na.rm = TRUE)
    nadd <- 5 * npos - nneg

    if (nadd >= 0) {
        vectadd <- rep(0, nadd)
        obsroc <- c(observedf, vectadd)
        indroc <- c(indexf, vectadd)
        ctrlval = 1

    } else {
        ctrlval = 0
    }
} else {
    obsroc <- observedf
    indroc <- indexf
    ctrlval = 1
}

# Writing averages of observation and model to output file
if (normalization == 2 && ctrlval == 1) {
    write(c(paste("AVGobs", as.character(round(mean(obsroc, na.rm = TRUE), 4)), sep = "\t"), 
        paste("AVGmodel", as.character(round(mean(indroc, na.rm = TRUE), 4)), sep = "\t")), 
        file = paste(prefix, "_results/", prefix, "_files/", prefix, "_averages.txt", 
            sep = ""), append = FALSE)
}

if (ctrlval == 1) {

    # Performing ROC analysis
    roc <- prediction(indexf, observedf)  #ROC prediction
    proc <- performance(roc, "tpr", "fpr")  #performance
    arocv <- performance(roc, "auc")  #area under curve

    # Formatting area under curve value for display
    arocb[i] = as.numeric(arocv@y.values)
    aroc[i] <- substitute(expression(iaa ~ aroca), list(iaa = format(ia), aroca = format(round(as.numeric(arocv@y.values), 
        digits = 3), nsmall = 3)))[2]

    # Plotting ROC curve
    plot(proc, add = TRUE, axes = F, xlab = NA, ylab = NA, lwd = lwidth[i], col = lcol[i], 
        lty = lltp[i])  #ROC curve
    mval = mval + 1

} else {
    arocb[i] <- NA
    aroc[i] <- NA
}

# Plotting bounding box and axes
box()
axis(side = 1, tck = -0.02, labels = NA)
axis(side = 2, tck = -0.02, labels = NA)
axis(side = 1, lwd = 0, line = -0.4)
axis(side = 2, lwd = 0, line = -0.4)

# Plotting axis labels
mtext(side = 1, "rFP/rON", line = 1.9)
mtext(side = 2, "rTP/rOP", line = 2)

# Plotting legend
aroch <- vector("expression", 1)
aroch[1] <- substitute(expression(AUROC))[2]
legend("bottomright", pch = NA, lty = c(0, lltp), lwd = c(0, lwidth), col = c("black", 
    lcol), legend = c(aroch[1], aroc), y.intersp = c(1, lintsp), text.col = c("black", 
    lcol), bty = "n", adj = c(0.1, 0.5))

# Writing area under curve to text file if ( ctrlval==1 ) { if ( normalization==1
# ) write(paste(prefix, '\t', round(as.numeric(arocv@y.values),3), sep=''),
# file='aucroc.txt', append=TRUE) if ( normalization==2 ) write(paste(prefix,
# '\t', round(as.numeric(arocv@y.values),3), sep=''), file='aucrocn.txt',
# append=TRUE) }

# Closing plot file
dev.off()
