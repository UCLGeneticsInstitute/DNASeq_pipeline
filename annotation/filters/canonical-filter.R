#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))

d <- as.data.frame(fread('file:///dev/stdin',showProgress=FALSE))

option_list <- list(
    make_option(c('--out'), default='', help='out file')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

d <- d[which(!is.na(d$SYMBOL)),]

d <- d[which(d$CANONICAL=='YES'),]

write.csv(d, file=opt$out, quote=FALSE, row.names=FALSE)


