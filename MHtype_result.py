#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: final_phasing_result.py
# @Author: willow
# @Site: 
# @Time: 6月 07, 2021
# ---
import argparse

import pandas as pd
from pathlib2 import Path
import os


class MHtype():
    def snp_sample(self, snpsample):
        df = pd.read_csv(snpsample, sep='\t')
        snpdict = dict(zip(df['ID'].tolist(), df['N'].tolist()))
        return snpdict

    def rawphasing(self, raw):
        df = pd.read_csv(raw, sep='\t')
        df['phasing'] = df['phasing'].str.replace("'", '')
        df['phasing'] = df['phasing'].str.replace("[", '')
        df['phasing'] = df['phasing'].str.replace("]", '')
        rawdict = dict(zip(df['2'].tolist(), df['phasing'].tolist()))
        return rawdict

    def mrginphsing(self, snpdict, phase):
        df1 = pd.read_csv(phase, sep='\t')
        df1['N'] = df1['id'].apply(lambda x: snpdict[x])
        df1['num'] = df1['N'].apply(lambda x: x.split('-')[0])
        phase1list = []
        list1 = df1.groupby(by='num')
        for k in list1.groups.keys():
            groupdf = list1.get_group(k)
            if len(groupdf) != 2:
                phase1list.append('none')
            else:
                if len(set(groupdf['basecount'].tolist()[0])) == 2 or len(set(groupdf['basecount'].tolist()[1])) == 2:
                    phase1list.append('none')
                    phase1list.append('none')
                else:
                    phase1list.append(
                        list(set(groupdf['basecount'].tolist()[0]))[0] + list(set(groupdf['basecount'].tolist()[1]))[0])
                    phase1list.append(
                        list(set(groupdf['basecount'].tolist()[0]))[0] + list(set(groupdf['basecount'].tolist()[1]))[0])
        return phase1list

    def phase1_phase2(self, phase1, phase2, snpdict):
        df1 = pd.read_csv(phase1, sep='\t')
        df2 = pd.read_csv(phase2, sep='\t')
        df1['phase1'] = self.mrginphsing(snpdict, phase1)
        df2['phase2'] = self.mrginphsing(snpdict, phase2)
        df = pd.merge(df1, df2, on='id')
        finaldf = pd.DataFrame()
        finaldf['id'] = df['id']
        finaldf['phase'] = df['phase1'] + ',' + df['phase2']
        finaldf['N'] = finaldf['id'].apply(lambda x: snpdict[x])
        finaldf['num'] = finaldf['N'].apply(lambda x: x.split('-')[0])
        return finaldf

    def finalphase(self, finaldf, rawdict, pileup):
        finaldf['tag'] = finaldf['id'] + '_' + finaldf['phase']
        finaldf['count'] = finaldf['num'].apply(lambda x: dict(finaldf['num'].value_counts())[x])
        finaldf = finaldf[finaldf['count'] == 2]
        finaldf['finalphase'] = self.findfinalphase(rawdict, finaldf['tag'].tolist(), pileup)
        del finaldf['count']
        del finaldf['phase']
        del finaldf['num']
        return finaldf

    def findfinalphase(self, rawdict, taglists, pileup):
        df = pd.read_csv(pileup, sep='\t')
        dictgenotype = dict(zip(df['id'].tolist(), df['basecount'].tolist()))
        varifyphase = []
        for i in range(0, len(taglists), 2):
            if 'none' in taglists[i]:
                varifyphase.append(rawdict[taglists[i].split('_')[0]])
                varifyphase.append(rawdict[taglists[i].split('_')[0]])
            else:
                varifyphase.append(taglists[i].split('_')[1])
                varifyphase.append(taglists[i].split('_')[1])
        return varifyphase

    def run(self, rawdir, raw, snpsamlpe, pileup):
        P = Path(rawdir)
        phase1 = P / 'H1/correct/phase1_correct_mpileup_results.txt'
        phase2 = P / 'H2/correct/phase2_correct_mpileup_results.txt'
        snpdict = self.snp_sample(snpsamlpe)
        rawdict = self.rawphasing(raw)
        df = self.phase1_phase2(phase1, phase2, snpdict)
        finaldf = self.finalphase(df, rawdict, pileup)
        finaldf.to_csv(P / 'finalphase.txt', sep='\t', index=False)
        os.system('rm -rf {}/H1'.format(P))
        os.system('rm -rf {}/H2'.format(P))

