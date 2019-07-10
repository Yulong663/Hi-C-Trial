# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 10:56:08 2017

@author: shaolab
"""

# Libraries
from matplotlib import pyplot as plt
from mirnylib.h5dict import h5dict
from matplotlib import cm
import os
import sys
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')


# Functions
def Heatmap(hdf5, sample, chrom, left_position, right_position):
    sfile = h5dict(hdf5, 'r')
    leftbound = int(left_position) / 40000
    rightbound = int(right_position) / 40000 + 1
    data = sfile[chrom + ' ' +
                 chrom][leftbound:rightbound, leftbound:rightbound]
    size = np.shape(data)[0]
    print size, np.shape(sfile[chrom + ' ' + chrom])
    #size = rightbound - leftbound
    plot_normalizer = matplotlib.colors.Normalize(vmin=1, vmax=10, clip=True)
    # plt.imshow(data,vmin=0,vmax=10,aspect='equal',cmap=cm.Reds)
    plt.imshow(
        data,
        norm=plot_normalizer,
        origin='lower',
        aspect='auto',
        cmap=cm.Reds)
    plt.gca().set_xlim(0, size)
    plt.gca().set_ylim(0, size)  # set the current axis length
    plt.colorbar(shrink=0.8)
    plt.gca().set_xticks(np.arange(0, size, 20))
    plt.gca().set_yticks(np.arange(0, size, 20))
    plt.gca().set_xticklabels(np.arange(leftbound, rightbound, 20))
    plt.gca().set_yticklabels(np.arange(leftbound, rightbound, 20))
    # plt.gca().set_xticklabels(np.arange(0,size,10))
    # plt.gca().set_xticklabels(np.arange(0,size,10))
    plt.setp(
        plt.gca().get_xticklabels(),
        rotation=90,
        ha="right",
        rotation_mode="anchor")
    plt.xlabel(chrom)
    plt.ylabel(str(left_position) + "-" + str(right_position))

    fig = plt.gcf()
    fig.set_size_inches(10, 8)
    fig.savefig(
        sample +
        "-" +
        str(left_position) +
        "-" +
        str(right_position) +
        ".png",
        dpi=300)
    # fig.clf()  clean all current axies but remain the current windows for
    # repeartly usage
    plt.close()

######
# NOTE: Before using the plot_real_eigen() function,check for the sign of eigenvalue !!!
######


def plot_real_eigen(chrom, left_pos, right_pos):

    # data process part
    file1 = open(chrom + ".txt", "r")
    eigenvector = np.array([])
    Chrms = [
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
    mono = [1, 1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1, 1, -1, -
            1, -1, 1, -1, 1, 1, -1, 1, -1, -1]  # mono_36 200K resolution
    # mac=[1, -1, 1, 1, -1, -1, 1, 1, -1, 1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1]          #mac_34 200K resolution
    # mono=[-1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, -1, -1, -1, 1, -1, -1, -1]        #THP-1 mono 200K resolution
    # mac=[1, -1, -1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, 1, -1, 1, -1, 1,
    # 1,-1, -1, 1, -1]         #THP-1 mac  200K resolution

    for line in file1:
        if line.strip().split()[0] == "NaN":
            #print line.strip().split()[0]
            eigenevector = np.append(eigenvector, 0)
        else:
            #print line.strip().split()[0]
            eigenvector = np.append(
                eigenvector, float(
                    line.strip().split()[0]))
    leftbound = left_pos / 200000
    rightbound = right_pos / 200000
    # check out for the specific list carefully!
    eigenvector = eigenvector[leftbound:rightbound] * mono[Chrms.index(chrom)]
    y1 = np.zeros(len(eigenvector))
    y2 = np.zeros(len(eigenvector))

    for i in xrange(len(eigenvector)):
        if eigenvector[i] >= 0:
            y1[i] = eigenvector[i]
        else:
            y2[i] = eigenvector[i]
    size = len(eigenvector)
    #temp=np.array(list(range(len(eigenvector)))) + leftbound
    temp = np.array(list(range(len(eigenvector))))

    # ploting part
    graph_name = chrom + "-" + \
        str(leftbound / 5.0) + "M-" + str(rightbound / 5.0) + "M"
    plt.plot(temp / 5.0, eigenvector, linewidth=0.1)
    plt.fill_between(temp / 5.0, y1, color="peru", alpha=1)
    plt.fill_between(temp / 5.0, y2, color="royalblue", alpha=1)
    # plt.xticks(np.arange(leftbound,rightbound,10))
    plt.xlabel('genome distance')
    plt.ylabel('Eigenvector value')

    plt.title(graph_name)
    # plt.savefig(graph_name+".png",dpi=300)
    plt.savefig(
        chrom +
        "-" +
        str(left_pos) +
        "-" +
        str(right_pos) +
        ".png",
        dpi=300)
    plt.close()


def Heat2Txt(hdf5, chrom):
    sfile = h5dict(hdf5, 'r')
    data = sfile[chrom + ' ' + chrom]
    np.savetxt(chrom + ".txt", data, delimiter="\t")


def all_heamap(hdf5):
    file1 = h5dict(hdf5, 'r')
    Chrms = [
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
    shape = {}
    for i in Chrms:
        for j in Chrms:
            key=i+' '+j
            data=file1[key]
            shape[key]=np.shape(data)
    size=0
    for name in shape.keys():
        size += shape[name][0]/24
    matrix = np.zeros((size,size),float)
    
# Main Function
if __name__ == '__main__':
    # Examples
    # hdf5='LiuXin_cHeatmap_Rep1_500k.hdf5'
    # chrom="chr1"
    # Heatmap(hdf5,chrom)
    # HeatmapCompare('FuDan_cHeatmap_500k.hdf5','LiuXin_cHeatmap_500k.hdf5','chr1')
    # Heat2Txt(hdf5,chrom)

    # Chrms=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    Chrms = [
        'chr2',
        'chr4',
        'chr3',
        'chr3',
        'chr4',
        'chr5',
        'chr6',
        'chr10',
        'chr6']
    left_position = [
        63979688,
        52791139,
        66000000,
        93600000,
        50000000,
        142000000,
        102000000,
        9000000,
        1000000]
    right_position = [
        87299687,
        76117138,
        90000000,
        115600000,
        58000000,
        153000000,
        114000000,
        21000000,
        11000000]
    #sample = sys.argv[2]
    #hdf5   = sys.argv[1]

    for i in xrange(len(Chrms)):
        # Heatmap(hdf5,sample,'chr2',63979688,87299687)
        # plot_eigen("chr2",63979688,87299687)
        # Heatmap(hdf5,sample,Chrms[i],left_position[i],right_position[i])
        plot_real_eigen(Chrms[i], left_position[i], right_position[i])
    # if not os.path.exists(sample):
        # os.mkdir(sample)
    #os.system('mv *png ' + sample + '/')
