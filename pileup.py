#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: new_pileup.py
# @Author: willow
# @Site: 
# @Time: 5月 06, 2021
# ---

import argparse

import pandas as pd
import re
import pysam


class Pileup():

    def LoadData(self, snprefile, pileupfile):
        df = pd.read_csv(snprefile, sep='\t', header=None)
        df['chr'] = df[0].apply(lambda x: 'chr{}'.format(int(x.split(".")[0][-2:])))
        df['new'] = df['chr'] + '_' + df[1].apply(lambda x: str(x))
        df1 = pd.read_csv(pileupfile, sep='\t', header=None)
        df1['new'] = df1[0] + '_' + df1[1].apply(lambda x: str(x))
        newdf = df1[df1['new'].isin(df['new'])]
        newdf['id'] = df[2].tolist()

        return newdf

    def SnpCall(self, newdf):
        newdf[4] = newdf[4].apply(lambda x: x.upper())
        newdf['basecount'] = self.Count(newdf[4].tolist(), newdf[3].tolist(), newdf[2].tolist())[0]
        newdf['maf'] = self.Count(newdf[4].tolist(), newdf[3].tolist(), newdf[2].tolist())[1]
        newdf['depth'] = self.Count(newdf[4].tolist(), newdf[3].tolist(), newdf[2].tolist())[2]
        return newdf

    def Linechange(self, line):
        line1 = list(line)
        for j in line1:
            if j == '+' or j == '-':
                i = line1.index(j)
                try:
                    if line1[i + 2] in ['A', 'T', 'C', 'G']:
                        nu = int(line1[i + 1])
                        del line1[i:i + nu + 2]
                    elif line1[i + 3] in ['A', 'T', 'C', 'G']:
                        nu = int(''.join(line1[i + 1:i + 3]))
                        del line1[i:i + nu + 2]

                except:
                    pass
        line1 = ''.join(line1)
        return line1

    def Count(self, line, length, rebase):
        bases = []
        maf = []
        depth = []
        for j in range(len(line)):
            # print(j)
            base = []
            pro = []
            base.append(rebase[j].upper())
            line1 = line[j]
            line1 = self.Linechange(line1)
            lengthnew = length[j] - line[j].count('*')
            # lengthnew = length[j]-line[j].count('*')-line[j].count('+')-line[j].count('-')
            # lengthnew = length[j]
            depth.append(lengthnew)
            pro.append(100 * (line1.count(".") + line1.count(",") + line1.count(rebase[j].upper())) / lengthnew)
            types = set(line1.upper())
            for i in types:
                if i in ["A", "T", "C", "G"] and i != rebase[j].upper():
                    base.append(i)
                    pro.append(100 * line1.count(i) / lengthnew)
            base_pro = dict(zip(base, pro))
            base_pro = sorted(base_pro.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
            try:
                maf.append(base_pro[1][1])
            except:
                maf.append([base_pro[0][1], len(base_pro)])
            if len(pro) != 1:
                try:
                    if base_pro[1][1] < 12:
                        bases.append(base_pro[0][0] * 2)
                    else:
                        bases.append(
                            base_pro[0][0] + base_pro[1][0])
                except:
                    bases.append(base_pro[0][0] * 2)

            else:
                bases.append(base_pro[0][0] * 2)

        return bases, maf, depth

    def run(self, snprefile, pileupfile, output):
        snp_pileup_df = self.LoadData(snprefile, pileupfile)
        SnpCallResult = self.SnpCall(snp_pileup_df)
        SnpCallResult.to_csv(output, index=False, sep='\t')
