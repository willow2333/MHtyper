#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: bamtofastaq.py
# @Author: willow
# @Site: 
# @Time: 5月 06, 2021
# ---
import argparse

import pandas as pd
import math
import os

class BamtoFastq():
    def clustering(self,file1,samplename):
        df1 = pd.read_csv(file1,sep='\t',header=None,quotechar="~")
        df1['tag'] = df1[3].apply(lambda x: math.floor(x/10000000))
        df1['new'] = df1[2]+'_'+df1['tag'].astype('str')
        df1 = df1[(df1[1] != 4) & (df1[9] != '*')]
        groups = df1.groupby(['new'])
        for i in range(len(list(groups))):
            newdf = '{}/{}.fastq'.format(samplename,i)
            newgroup = list(groups)[i][1]
            h = open(newdf, 'a')
            for j in range(len(newgroup)):
                h.write('@{}'.format(newgroup[0].tolist()[j]) + '\n')
                h.write(newgroup[9].tolist()[j] + '\n')
                h.write('+' + '\n')
                h.write(newgroup[10].astype('str').tolist()[j] + '\n')
            h.close()


