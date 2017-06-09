#! /bin/env python
from __future__ import print_function
import sys
import json
import argparse
import gzip
from pyfaidx import Fasta
import os
import pymongo
import json
import py2neo
from subprocess import call


parser=argparse.ArgumentParser(description='Arguments to annotation_importer.py')
parser.add_argument('--importdir', required=True)
parser.add_argument('--pheno4j_conn', required=True)
parser.add_argument('--pheno4j_user', required=True)
parser.add_argument('--pheno4j_password', required=True)
args=parser.parse_args()

basename='VEP'

Gene_filename='-'.join([basename,'Gene.csv'])
GeneticVariant_filename='-'.join([basename,'GeneticVariant.csv'])
gnomad_genomes_headers=['gnomad_genomes_AC_AFR', 'gnomad_genomes_AC_AMR', 'gnomad_genomes_AC_ASJ', 'gnomad_genomes_AC_raw', 'gnomad_genomes_AF_NFE', 'gnomad_genomes_AF_OTH', 'gnomad_genomes_AF_raw', 'gnomad_genomes_AN_AFR', 'gnomad_genomes_AC_EAS', 'gnomad_genomes_AC_Female', 'gnomad_genomes_AC_OTH', 'gnomad_genomes_AC_NFE', 'gnomad_genomes_AC_Male', 'gnomad_genomes_AC_FIN', 'gnomad_genomes_AF_AFR', 'gnomad_genomes_AF_AMR', 'gnomad_genomes_AF_ASJ', 'gnomad_genomes_AF_EAS', 'gnomad_genomes_AF_FIN', 'gnomad_genomes_AF_Female', 'gnomad_genomes_AF_Male', 'gnomad_genomes_AN_AMR', 'gnomad_genomes_AN_ASJ', 'gnomad_genomes_AN_EAS', 'gnomad_genomes_AN_FIN', 'gnomad_genomes_AN_Female', 'gnomad_genomes_AN_Male', 'gnomad_genomes_AN_NFE', 'gnomad_genomes_AN_OTH', 'gnomad_genomes_AN_raw', 'gnomad_genomes_Hom_AFR', 'gnomad_genomes_Hom_AMR', 'gnomad_genomes_Hom_ASJ', 'gnomad_genomes_Hom_EAS', 'gnomad_genomes_Hom_FIN', 'gnomad_genomes_Hom_Female', 'gnomad_genomes_Hom_Male', 'gnomad_genomes_Hom_NFE', 'gnomad_genomes_Hom_OTH', 'gnomad_genomes_Hom_raw', 'gnomad_genomes_Hom']
gnomad_exomes_headers=['gnomad_exomes_AC_AFR', 'gnomad_exomes_AC_AMR', 'gnomad_exomes_AC_ASJ', 'gnomad_exomes_AC_raw', 'gnomad_exomes_AF_NFE', 'gnomad_exomes_AF_OTH', 'gnomad_exomes_AF_raw', 'gnomad_exomes_AN_AFR', 'gnomad_exomes_AC_EAS', 'gnomad_exomes_AC_Female', 'gnomad_exomes_AC_OTH', 'gnomad_exomes_AC_NFE', 'gnomad_exomes_AC_Male', 'gnomad_exomes_AC_FIN', 'gnomad_exomes_AF_AFR', 'gnomad_exomes_AF_AMR', 'gnomad_exomes_AF_ASJ', 'gnomad_exomes_AF_EAS', 'gnomad_exomes_AF_FIN', 'gnomad_exomes_AF_Female', 'gnomad_exomes_AF_Male', 'gnomad_exomes_AN_AMR', 'gnomad_exomes_AN_ASJ', 'gnomad_exomes_AN_EAS', 'gnomad_exomes_AN_FIN', 'gnomad_exomes_AN_Female', 'gnomad_exomes_AN_Male', 'gnomad_exomes_AN_NFE', 'gnomad_exomes_AN_OTH', 'gnomad_exomes_AN_raw', 'gnomad_exomes_Hom_AFR', 'gnomad_exomes_Hom_AMR', 'gnomad_exomes_Hom_ASJ', 'gnomad_exomes_Hom_EAS', 'gnomad_exomes_Hom_FIN', 'gnomad_exomes_Hom_Female', 'gnomad_exomes_Hom_Male', 'gnomad_exomes_Hom_NFE', 'gnomad_exomes_Hom_OTH', 'gnomad_exomes_Hom_raw', 'gnomad_exomes_Hom']
GeneticVariant_headers=[ 'variantId', 'allele_string', 'start', 'end', 'seq_region_name', 'most_severe_consequence', 'strand', 'AC', 'allele_freq', 'AN', 'ExcessHet', 'FS', 'InbreedingCoeff', 'MLEAC', 'MLEAF', 'MQ', 'MQRankSum', 'ReadPosRankSum', 'VQSLOD', 'culprit', 'kaviar_AN', 'kaviar_AC', 'kaviar_AF', 'cadd_phred', 'cadd_raw']+gnomad_genomes_headers+gnomad_exomes_headers
TranscriptVariant_filename='-'.join([basename,'TranscriptVariant.csv'])
TranscriptVariant_headers=[ 'hgvsc','hgvsp','impact','exon','consequence_terms','fathmm_mkl_nc']
Transcript_filename='-'.join([basename,'Transcript.csv'])
# relationship files
GeneticVariantToGene_filename='-'.join([basename,'GeneticVariantToGene.csv'])
GeneticVariantToTranscriptVariant_filename='-'.join([basename,'GeneticVariantToTranscriptVariant.csv'])
TranscriptToGene_filename='-'.join([basename,'TranscriptToGene.csv'])

# sort the files and import into Pheno4J
py2neo.authenticate(args.pheno4j_conn, args.pheno4j_user, args.pheno4j_password)
graph = py2neo.Graph('http://%s/db/data/' % args.pheno4j_conn,secure=False,bolt=None, bolt_port=57687)


params=', '.join(['gv.%s=csvLine.%s' % (h,h) for h in GeneticVariant_headers])
s="""
LOAD CSV WITH HEADERS FROM "file:///%s" AS csvLine 
MERGE (gv:GeneticVariant {variantId: csvLine.variantId })
ON MATCH SET %s
ON CREATE SET %s
""" % (GeneticVariant_filename, params, params)
print(s)
#print(graph.run(s).stats())

s="""
LOAD CSV WITH HEADERS FROM "file:///%s" AS csvLine 
MERGE (tv:TranscriptVariant { hgvsc: csvLine.hgvsc, hgvsp: csvLine.hgvsp, impact: csvLine.impact, exon: csvLine.exon, consequence_terms: csvLine.consequence_terms, fathmm_mkl_nc: csvLine.fathmm_mkl_nc })
""" % TranscriptVariant_filename
print(s)
print(graph.run(s).stats())
# index

s="""
LOAD CSV WITH HEADERS FROM "file:///%s" AS csvLine 
MERGE (t:Transcript { transcript_id: csvLine.transcript_id })
""" % Transcript_filename
print(s)
print(graph.run(s).stats())

s="""
LOAD CSV WITH HEADERS FROM "file:///%s" AS csvLine 
MERGE (g:Gene { gene_id: csvLine.gene_id, gene_name:csvLine.gene_name })
""" % Gene_filename
print(s)
print(graph.run(s).stats())

# relationship

s="""
LOAD CSV WITH HEADERS FROM "file:///%s" AS csvLine 
MATCH (gv:GeneticVariant { variantId: csvLine.VariantId}),(p:Person { personId: csvLine.PersonId})
CREATE (gv)-[:HomVariantToPerson]->(person)
""" % hom_variant_filename
print(s)
print(graph.run(s).stats())

GeneticVariantToGene_filename
GeneticVariantToTranscriptVariant_filename
TranscriptToGene_filename





