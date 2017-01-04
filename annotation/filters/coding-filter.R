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
    make_option(c('--out'), help='out file')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

d <- d[which(!is.na(d$SYMBOL)),]

d <- d[which(d$CANONICAL=='YES'),]

d <- d[grep('transcript_ablation|splice_acceptor_variant|splice_donor_variant|stop_gained|frameshift_variant|stop_lost|start_lost|transcript_amplification|inframe_insertion|inframe_deletion|missense_variant|protein_altering_variant|splice_region_variant|incomplete_terminal_codon_variant|stop_retained_variant|synonymous_variant|coding_sequence_variant', d$Consequence),]

write.csv(d, file=opt$out, quote=FALSE, row.names=FALSE)


