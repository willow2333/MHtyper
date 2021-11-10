#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: run.py.py
# @Author: willow
# @Site: 
# @Time: 11月 02, 2021
# ---

import argparse
from bamtofastq import BamtoFastq
from MHtype_result import MHtype
from phasing import Phasing
from pileup import Pileup
import os


class Run():
    def FiltAligmnet(self, fastqfiles, reference, prefix):
        '''
            1.fastq files filtered by NanoFilt, the input files must be *.fq.gz
            2.alignment with minimap2
           '''
        commandFilter = 'gunzip -c {0} |NanoFilt  -q 7 -l 200 |gzip > {1}.fq.gz'
        commandAligment = 'minimap2 -t 10 -ax map-ont --secondary=no {0} {1}.fq.gz >{2}.sam'
        commandsam = 'samtools view -Sb {0}.sam >{1}.bam'
        commandsort = 'samtools sort -@6 -O bam -o {0}.sort.bam {1}.bam'
        commandindex = 'samtools index {0}.sort.bam'
        commandvcf = 'samtools mpileup -uvf {0} {1}.sort.bam | bcftools call -vm -Oz > bcftools.vcf.gz'
        os.system(commandFilter.format(fastqfiles, prefix))
        os.system(commandAligment.format(reference, prefix, prefix))
        os.system(commandsam.format(prefix, prefix))
        os.system(commandsort.format(prefix, prefix))
        os.system(commandindex.format(prefix))
        os.system(commandvcf.format(reference, prefix))

    def BedPileup(self, prefix):
        '''
            bed file & pileup file generation
        '''
        commandBed = 'bamToBed -i {0}.bam >{1}.bed'
        commandmpileup = 'samtools mpileup -f {0} -q 10 -Q 10 {1}.sort.bam >{2}_mpileup.txt'
        os.system(commandBed.format(prefix, prefix))
        os.system(commandmpileup.format(reference, prefix, prefix))

    def MHtype(self, truthvcf, prefix, marginpath, reference):
        '''
            MHtype resuts generation
        '''
        pileupdef = Pileup()
        phasingdef = Phasing()
        bamtofastqdef = BamtoFastq()
        MHtypedef = MHtype()
        pileupdef.run(truthvcf, '{}_mpileup.txt'.format(prefix), '{}_mpileup_results.txt'.format(prefix))
        phasingdef.run(truthvcf, '{}.sam'.format(prefix), '{}.bed'.format(prefix),
                       '{}_mpileup_results.txt'.format(prefix), '{}_phasing.txt'.format(prefix))
        commandmargin = '{0}/build/margin phase {1}.sort.bam {2} bcftools.vcf.gz {3}/params/misc/allParams.phase_vcf-default.json'
        os.system(commandmargin.format(marginpath, prefix, reference, marginpath))
        os.system('mkdir H1')
        os.system('mkdir H2')
        os.system('less {0}.sam |grep  ^[@].* >H1/phase1.bam'.format(prefix))
        os.system('samtools view output.haplotagged.bam | grep -E "HP:i:1"  >>H1/phase1.bam')
        os.system('less  {0}.sam |grep  ^[@].* >H2/phase2.bam'.format(prefix))
        os.system('samtools view output.haplotagged.bam | grep -E "HP:i:2"  >>H2/phase2.bam')
        os.chdir('./H1')
        os.system(
            "mkdir clustering;mkdir correction;samtools view phase1.bam |awk -v FS='\t' -v OFS='\t' '{if ($2 != 4 && $5 >= 10 ) print ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11)}'> phase1_new.txt")
        bamtofastqdef.clustering('phase1_new.txt', 'clustering')
        os.system('run_isoncorrect --t 15  --fastq_folder clustering  --outfolder correction')
        os.system('mkdir correct')
        os.chdir('./correct')
        os.system('cat ../correction/*/corrected_reads.fastq > phase1_correct.fastq')
        os.system(
            'minimap2 -t 10 -ax map-ont --secondary=no {0} phase1_correct.fastq >phase1_correct.sam'.format(reference))
        os.system('samtools view -Sb phase1_correct.sam >phase1_correct.bam')
        os.system('samtools sort -@6 -O bam -o phase1_correct.sort.bam phase1_correct.bam')
        os.system('samtools index phase1_correct.sort.bam')
        os.system('bamToBed -i phase1_correct.bam >phase1_correct.bed')
        os.system('samtools mpileup -f {0} -Q 10 phase1_correct.sort.bam >phase1_correct_mpileup.txt'.format(reference))
        pileupdef.run(truthvcf, 'phase1_correct_mpileup.txt', 'phase1_correct_mpileup_results.txt')
        os.chdir('../../')
        os.chdir('./H2')
        os.system(
            "mkdir clustering;mkdir correction;samtools view phase2.bam |awk -v FS='\t' -v OFS='\t' '{if ($2 != 4 && $5 >= 10 ) print ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11)}'> phase2_new.txt")
        bamtofastqdef.clustering('phase2_new.txt', 'clustering')
        os.system('run_isoncorrect --t 15  --fastq_folder clustering  --outfolder correction')
        os.system('mkdir correct')
        os.chdir('./correct')
        os.system('cat ../correction/*/corrected_reads.fastq > phase2_correct.fastq')
        os.system(
            'minimap2 -t 10 -ax map-ont --secondary=no {0} phase2_correct.fastq >phase2_correct.sam'.format(reference))
        os.system('samtools view -Sb phase2_correct.sam >phase2_correct.bam')
        os.system('samtools sort -@6 -O bam -o phase2_correct.sort.bam phase2_correct.bam')
        os.system('samtools index phase2_correct.sort.bam')
        os.system('bamToBed -i phase2_correct.bam >phase2_correct.bed')
        os.system('samtools mpileup -f {0} -Q 10 phase2_correct.sort.bam >phase2_correct_mpileup.txt'.format(reference))
        pileupdef.run(truthvcf, 'phase2_correct_mpileup.txt', 'phase2_correct_mpileup_results.txt')
        os.chdir('../../')
        MHtypedef.run(os.getcwd(), '{}_phasing.txt'.format(prefix), 'snp-sample.txt', '{}_mpileup_results.txt'.format(prefix))


if __name__ == '__main__':
    classrun = Run()
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastqfiles", help='The input *.fq.gz files.')
    parser.add_argument("--reference", help='The path of your ref.')
    parser.add_argument("--prefix", default='Test', help='The name of your Sample, default is "Test".')
    parser.add_argument("--truthvcf", help='The truth variant files in your research. ')
    parser.add_argument("--marginpath", help='The setup path of "Margin".')
    args = parser.parse_args()
    fastqfiles = args.fastqfiles
    reference = args.reference
    prefix = args.prefix
    truthvcf = args.truthvcf
    marginpath = args.marginpath
    classrun.FiltAligmnet(fastqfiles, reference, prefix)
    classrun.BedPileup(prefix)
    classrun.MHtype(truthvcf, prefix, marginpath, reference)
