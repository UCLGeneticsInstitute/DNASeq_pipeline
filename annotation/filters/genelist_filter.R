#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))

### Series of filters suggested by Adam.
message('*** coding FILTERING ***')

d <- as.data.frame(fread('file:///dev/stdin'))

option_list <- list(
    make_option(c('--gene.list'), help='list of genes to keep'),
    make_option(c('--out'), help='out file')
) 
option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

gene.list <- read.csv(opt$gene.list)[,1]

d <- d[which( d$SYMBOL %in% gene.list ),]

write.csv(d, file=opt$out, quote=FALSE, row.names=FALSE)

