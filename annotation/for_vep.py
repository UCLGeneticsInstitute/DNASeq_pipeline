#! /bin/env python
from __future__ import print_function
import sys

# zcat mainset_November2017.vcf.gz | cut -f-8 | python for_vep.py

def process_info(info,idx):
    info_string=[]
    for x in info.split(';'):
        if '=' not in x:
            info_string+=[x]
            continue
        k,v,=x.split('=')
        if ',' in v:
            info_string+=[k+'='+v.split(',')[idx]]
        else:
            info_string+=[k+'='+v]
    return ';'.join(info_string)

print('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO',sep='\t')
for l in sys.stdin:
    l=l.strip()
    if l.startswith('##'): continue
    if l.startswith('#CHROM'):
        headers=l.split('\t')
        continue
    d=l.split('\t')
    d=dict(zip(headers,d))
    id=d['ID']
    qual=d['QUAL']
    filter=d['FILTER']
    for idx,alt in enumerate(d['ALT'].split(',')):
        if alt=='*': alt='-'
        pos=int(d['POS'])
        ref=d['REF']
        info=process_info(d['INFO'],idx)
        print(d['#CHROM'],pos,id,ref,alt,qual,filter,info,sep='\t')
        if len(alt)>len(ref) and ref in alt :
            #insertion
            ind=alt.index(ref)
            pos=pos+ind+len(ref)-1
            alt=alt.replace(ref,'')
            ref='-'
            print(d['#CHROM'],pos,id,ref,alt,qual,filter,info,sep='\t')
        elif len(alt)<len(ref) and alt in ref :
            #deletion
            ind=ref.index(alt)
            pos=pos+ind+len(alt)
            ref=ref.replace(alt,'')
            alt='-'
            print(d['#CHROM'],pos,id,ref,alt,qual,filter,info,sep='\t')



