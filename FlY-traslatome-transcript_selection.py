import numpy as np
import re
from collections import OrderedDict
import sets
import sys


## INPUT  ：  the raw StringTie assemble transcript output
## OUTPUT:    the dict to store the gene and transcript information that satisfy the condition
## CONDITION: The gene have at least 2 transcript and fpkm of one transcript must > 2
def select_trans(input):
    file1=open(input,"r")
    transcript_dict=OrderedDict()
    transcript_dict1=OrderedDict()

    for lines in file1:
        line=lines.replace("\"","").replace(";","")
        line=re.split(r'[\s\t]',line)
        # print line
        gene=line[9]
        tr=line[11]
        cov=float(line[15])
        fpkm=float(line[17])
        tpm=float(line[19])
        tmp={'cov':cov,'fpkm':fpkm,'tpm':tpm}
        ## print tmp
        ## Use the 2 level nested dict structure to store the information
        if gene in transcript_dict:
            transcript_dict[gene][tr]=tmp
        else:
            transcript_dict[gene]={tr:tmp}
        #print len(transcript_dict)

    for gene in transcript_dict:
        if len(transcript_dict[gene])<2:
            continue
        else:
            for tr in transcript_dict[gene]:
                if transcript_dict[gene][tr]['fpkm']>2:
                    transcript_dict1[gene]=transcript_dict[gene]
                    break
                else:
                    continue

    # print len(transcript_dict),len(transcript_dict1),transcript_dict1[transcript_dict1.keys()[1]],transcript_dict1.keys()[1]                    
    return transcript_dict1



# begin to calculate the selection power of each transcript
def cal_RP(input_dict):
        tr_list=input_dict.keys()
        tr_fpkm=[]
        for key in tr_list:
            tr_fpkm.append(input_dict[key]['fpkm'])
        
        total_fpkm=float(sum(tr_fpkm))    
        tr_rp=[i/total_fpkm for i in tr_fpkm]
        for i in xrange(len(tr_list)):
            input_dict[tr_list[i]]['rp']=tr_rp[i]
        
        return input_dict
    

## remind that before calculate the Diff_CP for each pair of transcript between cyto and rib
## to sort it first and make sure each transcript in the pair is the same one
def cal_Diff_CP(cyto_trans_dict,rib_trans_dict,output):
    cyto=cal_RP(select_trans(cyto_trans_dict))
    rib=cal_RP(select_trans(rib_trans_dict))
    
    gene1=set(cyto.keys())
    gene2=set(rib.keys())
    
    file1=open(output,"w")
    common_gene=gene1.intersection(gene2)                          # figure out what is the common gene   
    
    ## output the gene ,transcript_ID,diff_cp,selection power
    for key in common_gene:
        diff_cp_list=[]
        tr_list1=set(cyto[key].keys())
        tr_list2=set(rib[key].keys())
        common_tr=tr_list1.intersection(tr_list2)                  # figure out what is the common gene   
        if len(common_tr) >1:                                      # To make sure that one gene has at least 2 transcript
            for tr in common_tr:
                diff_cp=cyto[gene][tr]['rp']-rib[gene][tr]['rp']   # diff_cp=cyto_cp -rib_cp
                diff_cp_list.append(diff_cp)
                positive_tran=min(diff_cp_list)                    # posotive selection means the minimal difference between cyto and rib is less
                negative_tran=max(diff_cp_list)                    # negative selection means the maximal difference
                selection_power=(abs(negative_tran) + abs(positive_tran))/2.0
                file1.write(gene + "\t" + tr + "\t" + str(diff_cp) + "\t" + str(diff_cp) + "\n")
                
        else:
            continue
            
            
if __name__ =="__main__":
    
    cal_Diff_CP(sys.argv[1],sys.argv[2],sys.argv[3])
            
