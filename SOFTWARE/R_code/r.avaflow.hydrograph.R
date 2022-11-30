#!/usr/bin/R

##############################################################################
#
# MODULE:       r.avaflow.hydrograph.R
# AUTHOR:       Martin Mergili
#
# PURPOSE:      The simulation model for avalanche and debris flows
#               Script for the creation of hydrograph plots
#
# COPYRIGHT:    (c) 2013 - 2021 by the author
#               (c) 2020 - 2021 by the University of Graz
#               (c) 2013 - 2020 by the BOKU University, Vienna
#               (c) 2015 - 2020 by the University of Vienna
#               (c) 1993 - 2021 by the R Development Core Team
#
# VERSION:      20210525 (25 May 2021)
#
#               This program is free software under the GNU General Public
#               License (>=v2).
#
##############################################################################

# Loading library
library(stats)

# Importing arguments defined in r.avaflow.py
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]  #prefix for file names
model <- as.integer(args[2])  #type of model
nhyd <- as.integer(args[3])  #number of hydrograph
nhydin <- as.integer(args[4])  #total number of input hydrographs
tint <- as.numeric(args[5])  #length of output time step

# Defining data file name
hydtable = paste(prefix, "_results/", prefix, "_files/", prefix, "_hydinfo", as.character(nhyd), 
    ".txt", sep = "")

# Creating vectors from data
thyd <- read.table(hydtable, skip = 1)[, 1]  #time passed
hhyd <- read.table(hydtable, skip = 1)[, 2]  #(PHASE 1) flow height
vhyd <- read.table(hydtable, skip = 1)[, 3]  #(PHASE 1) velocity
ehyd <- read.table(hydtable, skip = 1)[, 4]  #depth of (PHASE 1) entrainment/deposition
qhyd <- read.table(hydtable, skip = 1)[, 5]  #(PHASE 1) discharge

tmax <- round(max(thyd, na.rm = TRUE), 0)  #maximum value of time passed
hmax <- max(hhyd, na.rm = TRUE)  #maximum value of (PHASE 1) flow height
vmax <- max(abs(vhyd), na.rm = TRUE)  #maximum value of (PHASE 1) velocity
emax <- max(ehyd, na.rm = TRUE)  #maximum value of depth of (PHASE 1) entrainment/deposition
qmax <- max(abs(qhyd), na.rm = TRUE)  #maximum value of (PHASE 1) discharge

if (model > 3) {

    hhyd2 <- read.table(hydtable, skip = 1)[, 6]  #PHASE 2 flow height
    vhyd2 <- read.table(hydtable, skip = 1)[, 7]  #PHASE 2 velocity
    ehyd2 <- read.table(hydtable, skip = 1)[, 8]  #depth of PHASE 2 entrainment/deposition
    qhyd2 <- read.table(hydtable, skip = 1)[, 9]  #PHASE 2 discharge

    hmax2 <- max(hhyd2, na.rm = TRUE)  #maximum value of PHASE 2 flow height
    vmax2 <- max(abs(vhyd2), na.rm = TRUE)  #maximum value of PHASE 2 velocity
    emax2 <- max(ehyd2, na.rm = TRUE)  #maximum value of depth of PHASE 2 entrainment/deposition
    qmax2 <- max(abs(qhyd2), na.rm = TRUE)  #maximum value of PHASE 2 discharge

} else {
    hmax2 <- 0
    qmax2 <- 0
}

if (model == 7) {

    hhyd3 <- read.table(hydtable, skip = 1)[, 10]  #PHASE 3 flow height
    vhyd3 <- read.table(hydtable, skip = 1)[, 11]  #PHASE 3 velocity
    ehyd3 <- read.table(hydtable, skip = 1)[, 12]  #depth of PHASE 3 entrainment/deposition
    qhyd3 <- read.table(hydtable, skip = 1)[, 13]  #PHASE 3 discharge

    hmax3 <- max(hhyd3, na.rm = TRUE)  #maximum value of PHASE 3 flow height
    vmax3 <- max(abs(vhyd3), na.rm = TRUE)  #maximum value of PHASE 3 velocity
    emax3 <- max(ehyd3, na.rm = TRUE)  #maximum value of depth of PHASE 3 entrainment/deposition
    qmax3 <- max(abs(qhyd3), na.rm = TRUE)  #maximum value of PHASE 3 discharge

} else {
    hmax3 <- 0
    qmax3 <- 0
}

# Defining labels and scaling
ltext = vector("expression", 10)  #initializing label vector

if (nhyd <= nhydin) {
    textcol = "green"
    ltext[1] = substitute(expression(Input ~ hydrograph ~ IH ~ fnhyd), list(fnhyd = format(nhyd)))[2]  #title text for input hydrograph
} else {
    textcol = "purple"
    ltext[1] = substitute(expression(Output ~ hydrograph ~ OH ~ fnhyd), list(fnhyd = format(nhyd - 
        nhydin)))[2]  #title text for output hydrograph
}

ltext[2] = substitute(expression(Time ~ passed ~ italic(T) ~ (s)))[2]  #x axis label
ltext[3] = substitute(expression(Flow ~ height ~ italic(H) ~ (m)))[2]  #left y axis label

if (max(qmax, qmax2 + qmax3) < 1000) {
    # right y axis label:
    mconv = "1"
    ltext[4] = substitute(expression(Discharge ~ italic(Q) ~ (m^3 ~ s^-1)))[2]
} else if (max(qmax, qmax2 + qmax3) < 1e+06) {
    mconv = "0.001"
    ltext[4] = substitute(expression(Discharge ~ italic(Q) ~ (10^3 ~ m^3 ~ s^-1)))[2]
} else {
    mconv = "0.000001"
    ltext[4] = substitute(expression(Discharge ~ italic(Q) ~ (10^6 ~ m^3 ~ s^-1)))[2]
}
qmax <- qmax * as.numeric(mconv)  #scaling of right y axis
qmax2 <- qmax2 * as.numeric(mconv)
qmax3 <- qmax3 * as.numeric(mconv)

# Creating plot file
hydplot = paste(prefix, "_results/", prefix, "_plots/", prefix, "_hydrograph", nhyd, 
    ".png", sep = "")
png(filename = hydplot, width = 15, height = 10, units = "cm", res = 300)

# Defining external margins
par(mar = c(3.2, 3.2, 1, 3.2))

# Defining internal margins
if (((hmax + hmax2) > (2 * hmax3)) && ((qmax + qmax2) > (2 * qmax3))) {
    marg <- 1.1
    tpos <- -(hmax + hmax2) * 0.85
    lpos <- "bottomright"
} else if ((hmax2 + hmax3) > 2 * hmax && (qmax2 + qmax3) > 2 * qmax) {
    marg <- 1.1
    tpos <- hmax3 * 0.85
    lpos <- "topright"
} else {
    marg <- 1.3
    lpos <- "bottomright"
    tpos <- max(hmax + hmax2, hmax3) * 1.25
}

# Plotting bars
if (model <= 3) {
    barplot(qhyd * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(-qmax * 
        marg, qmax * marg), border = NA, col = rgb(0.5, 0.85, 0.5))  #discharge
} else if (model > 3 && model < 7) {
    barplot(qhyd * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(-max(qmax2 * 
        marg, qmax * marg), max(qmax2 * marg, qmax * marg)), border = NA, col = rgb(0.45, 
        0.3, 0, 0.5))  #PHASE 1 discharge
} else if (model == 7) {
    barplot(qhyd * as.numeric(mconv) + qhyd2 * as.numeric(mconv), axes = F, xlab = NA, 
        ylab = NA, ylim = c(-max(qmax3 * marg, qmax2 * marg + qmax * marg), max(qmax3 * 
            marg, qmax2 * marg + qmax * marg)), border = NA, col = rgb(0.6, 0, 0, 
            0.5))  #PHASE 1 discharge
}

par(new = TRUE)
if (model == 7) {
    barplot(-qhyd3 * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(-max(qmax3 * 
        marg, qmax2 * marg + qmax * marg), max(qmax3 * marg, qmax2 * marg + qmax * 
        marg)), border = NA, col = rgb(0, 0, 0.6, 0.5))  #PHASE 3 discharge
}

par(new = TRUE)
if (model > 3 && model < 7) {
    barplot(-qhyd2 * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(-max(qmax2 * 
        marg, qmax * marg), max(qmax2 * marg, qmax * marg)), border = NA, col = rgb(0, 
        0, 0.6, 0.5))  #PHASE 2 discharge
} else if (model == 7) {
    barplot(qhyd2 * as.numeric(mconv), axes = F, xlab = NA, ylab = NA, ylim = c(-max(qmax3 * 
        marg, qmax2 * marg + qmax * marg), max(qmax3 * marg, qmax2 * marg + qmax * 
        marg)), border = NA, col = rgb(0, 0.6, 0, 0.5))  #PHASE 2 discharge
}

axis(side = 4, tck = -0.02, labels = NA)  #y axis for velocity and random kinetic energy
axis(side = 4, lwd = 0, line = -0.4)

# Plotting lines
par(new = TRUE)
if (model <= 3) {
    plot(x = thyd, y = hhyd, lty = 1, lwd = 2, col = rgb(0.45, 0.3, 0), axes = F, 
        xlab = NA, ylab = NA, type = "l", ylim = c(-hmax * marg, hmax * marg))  #flow height
} else if (model > 3 && model < 7) {
    plot(x = thyd, y = hhyd, lty = 1, lwd = 2, col = rgb(0.45, 0.3, 0), axes = F, 
        xlab = NA, ylab = NA, type = "l", ylim = c(-max(hmax * marg, hmax2 * marg, 
            hmax3 * marg), max(hmax * marg, hmax2 * marg, hmax3 * marg)), xlim = c(-tint/2, 
            tmax + tint/2))  #PHASE 1 flow height
    lines(x = thyd, y = -hhyd2, lty = 1, lwd = 2, col = rgb(0, 0, 0.8), xlim = c(-tint/2, 
        tmax + tint/2))  #PHASE 2 flow height
} else if (model == 7) {
    plot(x = thyd, y = hhyd, lty = 1, lwd = 2, col = rgb(0.6, 0, 0), axes = F, xlab = NA, 
        ylab = NA, type = "l", ylim = c(-max(hmax * marg, hmax2 * marg, hmax3 * marg), 
            max(hmax * marg, hmax2 * marg, hmax3 * marg)), xlim = c(-tint/2, tmax + 
            tint/2))  #PHASE 1 flow height
    lines(x = thyd, y = hhyd2, lty = 1, lwd = 2, col = rgb(0, 0.6, 0), xlim = c(-tint/2, 
        tmax + tint/2))  #PHASE 2 flow height
}

if (model == 7) {
    lines(x = thyd, y = -hhyd3, lty = 1, lwd = 2, col = rgb(0, 0, 0.6), xlim = c(-tint/2, 
        tmax + tint/2))  #PHASE 3 flow height
}
abline(0, 0, col = "black", lwd = 1)  #zero line

# Text and legends
text(x = tmax/2, y = tpos, labels = ltext[1], font = 2, col = textcol, cex = 1.1)  #title text

if (model <= 3) {

    ltext[5] = substitute(expression(italic(H)))[2]  #flow height
    ltext[6] = substitute(expression(italic(Q)))[2]  #discharge

    legend(lpos, legend = ltext[5], lty = 1, lwd = 2, col = "brown", text.col = "black", 
        bty = "n", inset = c(1.025, -0.025), horiz = TRUE)
    # flow height
    legend(lpos, legend = ltext[6], border = NA, fill = rgb(0.5, 0.85, 0.5), text.col = "black", 
        bty = "n", inset = c(0.225, -0.025), horiz = TRUE)
    # discharge

} else if (model > 3 && model < 7) {

    ltext[5] = substitute(expression(italic(H)[s]))[2]  #PHASE 1 flow height
    ltext[6] = substitute(expression(italic(Q)[s]))[2]  #PHASE 1 discharge
    ltext[7] = substitute(expression(italic(-H)[f]))[2]  #PHASE 2 flow height
    ltext[8] = substitute(expression(italic(-Q)[f]))[2]  #PHASE 2 discharge

    legend(lpos, legend = ltext[5], lty = 1, lwd = 2, col = rgb(0.45, 0.3, 0), text.col = "black", 
        bty = "n", inset = c(0.825, -0.025), horiz = TRUE)
    # PHASE 1 flow height
    legend(lpos, legend = ltext[7], lty = 1, lwd = 2, col = rgb(0, 0, 0.8), text.col = "black", 
        bty = "n", inset = c(0.625, -0.025), horiz = TRUE)
    # PHASE 2 flow height
    legend(lpos, legend = ltext[6], border = NA, fill = rgb(0.45, 0.3, 0, 0.5), text.col = "black", 
        bty = "n", inset = c(0.225, -0.025), horiz = TRUE)
    # PHASE 1 discharge
    legend(lpos, legend = ltext[8], border = NA, fill = rgb(0, 0, 0.6, 0.5), text.col = "black", 
        bty = "n", inset = c(0.025, -0.025), horiz = TRUE)
    # PHASE 2 discharge
} else if (model == 7) {

    par(cex = 0.75)

    ltext[5] = substitute(expression(italic(H)[P1]))[2]  #PHASE 1 flow height
    ltext[6] = substitute(expression(italic(Q)[P1]))[2]  #PHASE 1 discharge
    ltext[7] = substitute(expression(italic(H)[P2]))[2]  #PHASE 2 flow height
    ltext[8] = substitute(expression(italic(Q)[P2]))[2]  #PHASE 2 discharge
    ltext[9] = substitute(expression(italic(-H)[P3]))[2]  #PHASE 3 flow height
    ltext[10] = substitute(expression(italic(-Q)[P3]))[2]  #PHASE 3 discharge

    legend(lpos, legend = ltext[5], lty = 1, lwd = 2, col = rgb(0.6, 0, 0), text.col = "black", 
        bty = "n", inset = c(0.85, -0.025), horiz = TRUE)
    # PHASE 1 flow height
    legend(lpos, legend = ltext[7], lty = 1, lwd = 2, col = rgb(0, 0.6, 0), text.col = "black", 
        bty = "n", inset = c(0.7, -0.025), horiz = TRUE)
    # PHASE 2 flow height
    legend(lpos, legend = ltext[9], lty = 1, lwd = 2, col = rgb(0, 0, 0.6), text.col = "black", 
        bty = "n", inset = c(0.55, -0.025), horiz = TRUE)
    # PHASE 3 flow height
    legend(lpos, legend = ltext[6], border = NA, fill = rgb(0.6, 0, 0, 0.5), text.col = "black", 
        bty = "n", inset = c(0.35, -0.025), horiz = TRUE)
    # PHASE 1 discharge
    legend(lpos, legend = ltext[8], border = NA, fill = rgb(0, 0.6, 0, 0.5), text.col = "black", 
        bty = "n", inset = c(0.2, -0.025), horiz = TRUE)
    # PHASE 2 discharge
    legend(lpos, legend = ltext[10], border = NA, fill = rgb(0, 0, 0.6, 0.5), text.col = "black", 
        bty = "n", inset = c(0.05, -0.025), horiz = TRUE)
    # PHASE 3 discharge
    par(cex = 1)

}

# Plotting bounding box and axes
box()
axis(side = 1, tck = -0.02, labels = NA, at = c(0, tmax/10, 2 * tmax/10, 3 * tmax/10, 
    4 * tmax/10, 5 * tmax/10, 6 * tmax/10, 7 * tmax/10, 8 * tmax/10, 9 * tmax/10, 
    tmax))
# x axis
axis(side = 1, lwd = 0, line = -0.4, at = c(0, tmax/10, 2 * tmax/10, 3 * tmax/10, 
    4 * tmax/10, 5 * tmax/10, 6 * tmax/10, 7 * tmax/10, 8 * tmax/10, 9 * tmax/10, 
    tmax))
axis(side = 2, tck = -0.02, labels = NA)  #y axis for lines
axis(side = 2, lwd = 0, line = -0.4)

# Plotting axis labels
mtext(side = 1, ltext[2], line = 2)  #x axis
mtext(side = 2, ltext[3], line = 2)  #y axis for lines
if (model <= 3) {
    mtext(side = 4, ltext[4], line = 2)  #y axis for bars
} else if (model > 3 && model < 7) {
    mtext(side = 4, ltext[4], line = 2)  #y axis for bars
} else if (model == 7) {
    mtext(side = 4, ltext[4], line = 2)  #y axis for bars
}

# Closing plot file
dev.off()
