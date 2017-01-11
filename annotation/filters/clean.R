#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))

d <- as.data.frame(fread('file:///dev/stdin',showProgress=FALSE))

option_list <- list(
    make_option(c('--headers'), default='headers.txt', help='two column header file, only keep headers defined in this file.'),
    make_option(c('--out'), default='', help='out file')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

message('*** CLEAN ***')

#v <- v[,-which(colnames(v) %in% c('DISTANCE','CADD_RAW','CADD_PHRED','SOMATIC','CLIN_SIG','Feature_type'))]
#d <- d[,-which(colnames(d) %in% c('SYMBOL_SOURCE','HGNC_ID','CANONICAL','CADD_PHRED','GO','MISS','WT','AF','SOMATIC','ID'))]
#v <- v[,-grep('ESP',colnames(v))]
#v <- v[,which(!grepl('SYMBOL_SOURCE|HGNC_ID|CANONICAL|CADD_PHRED|GO|MISS|WT|AF',colnames(v)))]
#v <- v[,which(!grepl('geno',colnames(v)))]
#v <- v[,which(!grepl('depth',colnames(v)))]
#
#cleanup GO terms (is done in the GO filter)

# move HET HOM to front
# move kaviar to front

if (file.exists(opt$headers)) {
    headers <- as.data.frame(fread(opt$headers))
}

#cols <- intersect(colnames(d),headers$name)
cols <- headers$name

message(setdiff(headers$name,colnames(d)))

write.csv(d[,cols], file=opt$out, quote=FALSE, row.names=FALSE)
