#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: phasing.py
# @Author: willow
# @Site: 
# @Time: 3月 25, 2021
# ---
import argparse

import pandas as pd
import csv
import re
import pysam


class Phasing():

    def LoadData(self, pilefile, snprefile, bedfile, samfile):
        piledf = pd.read_csv(pilefile, sep='\t', header=None)
        piledf['new'] = piledf[0] + '_' + piledf[1].astype('str')
        new_genotype = dict(zip(piledf['new'].tolist(), piledf[8].tolist()))

        snpdf = pd.read_csv(snprefile, sep='\t', header=None)
        snpdf['chr'] = snpdf[0].apply(lambda x: 'chr{}'.format(int(x.split(".")[0][-2:])))
        snpdf['new'] = snpdf['chr'] + '_' + snpdf[1].apply(lambda x: str(x))
        # snpdf = snpdf[snpdf['new'].isin(piledf['new'])]
        snpdf['type'] = snpdf['new'].apply(lambda x: new_genotype[x])

        beddf = pd.read_csv(bedfile, sep='\t', header=None)

        samdf = pd.read_csv(samfile, sep='\t', header=None, comment='@', names=[i for i in range(25)],
                            quoting=csv.QUOTE_NONE)
        samdf = samdf[samdf[9] != '*']
        samdf = samdf[samdf[4] > 0]
        samdf['chr_pos_id'] = samdf[2] + '_' + samdf[3].astype('str') + '_' + samdf[0]

        return snpdf, beddf, samdf

    def Newseq(self, samdf):
        seqs = samdf[9].tolist()
        cigars = samdf[5].tolist()
        finalseqs = []
        for i in range(len(seqs)):
            finalseq = ''
            seq = seqs[i]
            citgar = cigars[i]
            citgar_list = []
            nu = 0
            for j in range(len(citgar)):
                if citgar[j] not in ['M', 'D', 'S', 'H', 'I']:
                    continue
                else:
                    citgar_list.append(citgar[nu:j + 1])
                    nu = j + 1
            nu1 = 0
            for h in citgar_list:
                if h[-1] == 'M':
                    # if h[-1] not in ['I','D','H','S']:
                    finalseq += seq[nu1:int(h[0:-1]) + nu1]
                    nu1 = nu1 + int(h[0:-1])
                elif h[-1] == 'D':
                    finalseq += '*' * (int(h[0:-1]))
                elif h[-1] == 'S':
                    nu1 = nu1 + int(h[0:-1])
                elif h[-1] == 'I':
                    nu1 = nu1 + int(h[0:-1])
            finalseqs.append(finalseq)

        samdf['new_seq'] = finalseqs
        return samdf

    def Convert(self, snpdf, beddf, samdf):
        new = snpdf[1].tolist()
        chrs = snpdf['new'].tolist()
        snptypes = snpdf['type'].tolist()
        finaltypeslist = []
        #finalreads = []
        for i in range(0, len(new), 2):
            totaltype1 = snptypes[i]
            totaltype2 = snptypes[i + 1]
            start = new[i]
            end = new[i + 1]
            newdf1 = beddf[(beddf[1].astype('int') < int(start)) & (int(end) < beddf[2].astype('int'))]
            newdf1['1_new'] = newdf1[1].apply(lambda x: str(int(x) + 1))
            newdf1['chr_pos_id'] = newdf1[0] + '_' + newdf1['1_new'] + '_' + newdf1[3]
            finaldf = samdf[(samdf['chr_pos_id'].isin(newdf1['chr_pos_id'])) & (samdf[2] == chrs[i].split('_')[0])]
            finaldf = self.Newseq(finaldf)
            types = []
            finaltypes = {}
            rawtypes = {}
            pos = finaldf[3].tolist()
            seqs = finaldf['new_seq'].tolist()
            for j in range(len(seqs)):
                snp1 = int(start) - pos[j]
                snp2 = int(end) - pos[j]
                try:
                    types.append(seqs[j][snp1] + seqs[j][snp2])
                except:
                    # print(start)
                    pass
            for k in set(types):
                if k[0] in totaltype1 and k[1] in totaltype2:
                    finaltypes[k] = types.count(k)
                    rawtypes[k] = types.count(k)
                else:
                    rawtypes[k] = types.count(k)

            a = sorted(finaltypes.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)

            if len(a) > 1:
                k1 = totaltype1.replace(a[0][0][0], '') if len(set(totaltype1)) == 2 else totaltype1[0]
                k2 = totaltype2.replace(a[0][0][1], '') if len(set(totaltype2)) == 2 else totaltype2[0]
                type2 = k1 + k2
                finaltypeslist.append([a[0][0], type2])
                finaltypeslist.append([a[0][0], type2])

            elif len(a) == 1:
                finaltypeslist.append(a[0][0])
                finaltypeslist.append(a[0][0])

            else:
                finaltypeslist.append('')
                finaltypeslist.append('')
        snpdf['phasing'] = finaltypeslist
        #snpdf['reads'] = finalreads

        return snpdf

    def run(self, snprefile, samfile, bedfile, pilefile, output):
        snpdf, beddf, samdf = self.LoadData(pilefile, snprefile, bedfile, samfile)
        finaldf = self.Convert(snpdf, beddf, samdf)
        finaldf.to_csv(output, index=False, sep='\t')
