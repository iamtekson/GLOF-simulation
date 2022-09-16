#!/usr/bin/R

##############################################################################
#
# MODULE:       r.avaflow.merge.R
# AUTHOR:       Martin Mergili
#
# PURPOSE:      Script for merging of parameter and
#               evaluation tables
#
# COPYRIGHT:    (c) 2017 - 2021 by the author
#               (c) 2020 - 2021 by the University of Graz
#               (c) 2017 - 2020 by the BOKU University, Vienna
#               (c) 2017 - 2020 by the University of Vienna
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

# Importing arguments (defined in r.ranger.py)
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]

# Defining input parameter file name
inparam = paste(prefix, "_results/", prefix, "_files/", prefix, "_params.txt", sep = "")

# Defining evaluation parameter file name
invalid = paste(prefix, "_results/", prefix, "_files/", prefix, "_evaluation.txt", 
    sep = "")

parpar <- read.table(inparam, header = TRUE)
parval <- read.table(invalid, header = TRUE)
parcomb <- merge(parpar, parval, by.x = 1, by.y = 1)  #combining tables
write.table(parcomb, file = paste(prefix, "_results/", prefix, "_files/", prefix, 
    "_allparams0.txt", sep = ""), append = FALSE, quote = FALSE, sep = "\t")
