# -*- coding:utf8 -*-

import os
import pymol
import numpy as np
import sys

def hairpin(hicfile,chrom,window1=2,window2=20):
    if type(hicfile)!=str:
        data=np.asarray(hicfile)
    else:
        if hicfile.endswith('txt'):
            data=np.genfromtxt(hicfile)
        else:
            if hicfile.endswith('npy'):
                data=np.load(hicfile)
            else:
                raise IOError('Wrong input format!')
    N=len(data)
    signals=np.zeros(N)
    for i in range(N):
        if i<window2/2 or i>=N-window2/2:
            continue
        matrix=data[i-window2/2:i+window2/2+1,i-window2/2:i+window2/2+1]
        tmp=np.arange(len(matrix))
        index1=(tmp[:,None]+tmp[None,:]>=((window2/2-window1/2)*2))*\
            (tmp[:,None]+tmp[None,:]<=((window2/2+window1/2)*2))
        index2=(tmp[:,None]+tmp[None,:]<((window2/2-window1/2)*2))
        index3=(tmp[:,None]+tmp[None,:]>((window2/2+window1/2)*2))
        ridge=np.mean(matrix[index1])
        left=np.mean(matrix[index2])
        right=np.mean(matrix[index3])
        if left==0 or right==0:
            continue
        signal=ridge/np.sqrt(left*right)
        signals[i]=signal
    # return signals
    position2bed(signals,chrom)

def position2bed(vector, chrm):
    for i in xrange(len(vector)):
        print chrm + "\t" + str(i*40000) +str(((i+1)*40000)-1) + "\t" + str(vector[i])

def find_hairpin_region(vector,cutoff=0.95):
    index1=vector >= cutoff
    hairpin_region1=np.array([])
    for i in xrange(len(vector)):
        if inde1[i] == True:
            hairpin_region=np.append(hairpin_region1,i)
        else:
            continue

    hairpin_region2={}
    count=1
    for j in xrange(len(hairpin_region1)):
        if hairpin_region1[j]<= hairpin_region1[j+1] :
            if count in hairpin_region2:
                hairpin_region2[count]=np.append(hairpin_region2[count],hairpin_region1[j])
            else:
                hairpin_region2[count]=np.array(hairpin_region2[j])
        elif hairpin_region1[j] > hairpin_region1[j+1] and hairpin_region1[j]+1 == hairpin_region1[j+1]:
            hairpin_region2[count]=np.append(hairpin_region2[count],hairpin_region1[j])

        else:
            count += 1
    local_optimal=np.array([])
    for key in hairpin_region2.keys():


if __name__ == '__main__':
    Chrms=[
        'chr1',
        'chr2',
        'chr3',
        'chr4',
        'chr5',
        'chr6',
        'chr7',
        'chr8',
        'chr9',
        'chr10',
        'chr11',
        'chr12',
        'chr13',
        'chr14',
        'chr15',
        'chr16',
        'chr17',
        'chr18',
        'chr19',
        'chr20',
        'chr21',
        'chr22',
        'chrX',
	    'chrY']

    #dict1={}
    for i in Chrms:
        hairpin(i + '.txt', window1=2, window2=20,chrom=i)
        #dict1[i]=hairpin(i+'.txt',window1=2,window2=20)
    #for j in dict1:
        #np.savetxt(j+'_hairpin.txt',dict1[j],delimiter='\t',fmt='%f')





