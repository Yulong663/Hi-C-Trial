# -*- coding:utf8 -*-

import os
import pymol
import numpy as np
import sys

def hairpin(hicfile,window1=2,window2=20,chrom):
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
    #position2bed(signals,chrom)
    find_hairpin_region(signals,cutoff=0.95,chrom)

def position2bed(vector, chrm):
    for i in xrange(len(vector)):
        print chrm + "\t" + str(i*40000) +str(((i+1)*40000)-1) + "\t" + str(vector[i])

def find_hairpin_region(vector,chrom,cutoff=0.95):
    hairpin_pos=np.array([])
    vector = np.array(vector)
    for i in xrange(len(vector)):
        if i == len(vector) -1 or i ==0:
            continue
        hairpin_pos=np.append(hairpin_pos,i)

    ### Every pair of key:value in hairpin_region2 was used to store a piece of region which signal value above the cutoff
    ### The hairpin_region1 was used to store the position information of regions that signal value >= cutoff
    hairpin_region={}
    count=1
    hairpin_region[count]=np.array([])
    for j in hairpin_pos:
        if vector[j]<= vector[j+1] :
            if vector[j-1] <= vector[j]:
                hairpin_region[count] = np.append(hairpin_region[count],j)
            else:
                count +=1
                hairpin_region[count]=np.array([j])
        else :
            hairpin_region[count]=np.append(hairpin_region[count],[j])

    ### to remove those regions that maximum of the peak less than cutoff (i.e remove the noise)
    for key in hairpin_region.keys():
        if np.max(hairpin_region[key]) <= cutoff:
            continue
        else：
            left_bound=hairpin_region[key][0]
            right_bound=hairpin_region[key][-1]
            print chrom + "\t" + str(int(left_bound)*40000) + "\t" +str((int(right_bound)*40000)-1)

		
		
		

def find_hairpin_region_1(vector,chrom,cutoff=0.85):
    vector=np.array(vector)
    hairpin_region1=np.array([])
    for i in xrange(len(vector)):
        hairpin_region1=np.append(hairpin_region1,i)                # store the position information of hairpin_signal vector

    hairpin_region2={}
    count=1
    for j in xrange(len(hairpin_region1)):
        if j == len(hairpin_region1) -1:
            continue

        if vector[hairpin_region1[j]] <= vector[hairpin_region1[j+1]] and hairpin_region1[j]+1 == hairpin_region1[j+1]:
            if count in hairpin_region2:
                hairpin_region2[count]=np.append(hairpin_region2[count],hairpin_region1[j])
            else:
                hairpin_region2[count]=np.array([hairpin_region1[j]])
        ### under this condition test,every piece of hairpin_region will lost the final block.If not to do so , the pieces will connect unwillingly !
        elif vector[hairpin_region1[j]] > vector[hairpin_region1[j+1]] and hairpin_region1[j]+ 1 == hairpin_region1[j+1]:
            if count in hairpin_region2:
                hairpin_region2[count]=np.append(hairpin_region2[count],hairpin_region1[j])
            else:
                hairpin_region2[count]=np.array([hairpin_region1[j]])                ### remind if there is only 1 point, then the region will a np.array(value), and report wrong message
        else:
            count += 1
    local_optimal=np.array([])
    for key in hairpin_region2.keys():
        local_optimal=np.append(local_optimal,np.max(hairpin_region2[key]))
        left_bound=hairpin_region2[key][0]
        right_bound=hairpin_region2[key][len(hairpin_region2[key])-1] + 1
        print chrom + "\t" + str(int(left_bound)*40000) + "\t" +str((int(right_bound) +1) * 40000 - 1)


def old_find_hairpin_region(vector,cutoff=0.95,chrom):
    hairpin_region1=np.array([])
    index1= vector >= cutoff
    for i in xrange(len(vector)):
        if i == len(vector) -1:
            continue

        if index1[i] == True:
            hairpin_region1=np.append(hairpin_region1,i)
        else:
            continue
    ### Every pair of key:value in hairpin_region2 was used to store a piece of region which signal value above the cutoff
    ### The hairpin_region1 was used to store the position information of regions that signal value >= cutoff
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
        local_optimal=np.append(local_optimal,np.max(hairpin_region2[key]))
        left_bound=hairpin_region2[key][0]
        right_bound=hairpin_region2[key][len(hairpin_region2[key]) - 1]
        print chrom + "\t" + str(int(left_bound)*40000) + "\t" +str((int(right_bound)*40000)-1)


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
        hairpin(i + '.txt', chrom=i,window1=2, window2=20)
        #dict1[i]=hairpin(i+'.txt',window1=2,window2=20)
    #for j in dict1:
        #np.savetxt(j+'_hairpin.txt',dict1[j],delimiter='\t',fmt='%f')





