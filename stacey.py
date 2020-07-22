# -*- coding: utf-8 -*-
"""
Read crystallines large file
"""

#%%
#Importing libraries...
from Bio.SeqIO import parse 
#from Bio.SeqRecord import SeqRecord 
#from Bio.Seq import Seq 
from Bio.SeqUtils.ProtParam import ProteinAnalysis

import scipy as sc
import pylab as pl
import matplotlib.pyplot as plt
import numpy as np
#%% find all elements in alist that match keyw, case INsensitive 
def makemask(alist,keyw):
    mask=np.zeros(len(alist),bool)
    for i in range(len(alist)):
        mask[i]= ( alist[i].upper().find(keyw.upper()) >=0 )
    return mask
        
#%%
file=open("/home/mysh/colours/scripts/Crystallins/txt/allgenbank.txt")
#%%
#Creating the list of dictionaries with percentages
sequence_list=[]
names_list=[]
sources_list = []
desc_list = []
taxo_list = []
keyw_list = []
for record in  parse(file, "genbank"):
   cdsnum=0
   for feat in record.features:
       if(feat.type == 'CDS'): #CDS is one of the features (conserved domain sequence)
           try: #throws an error if there is a stop codon in the nucleotide sequence
               prot=feat.translate(record.seq)
               analysed_seq = ProteinAnalysis(str(prot)) #creating another fucking class ProteinAnalysis
               sequence_list.append(analysed_seq.get_amino_acids_percent()) #invoking method on this class, it returns a dictionary, we store it in the list
               #saving lists for the other information we might need
               names_list.append(str(record.name)+ "_CDS#" + str(cdsnum))                    
               sources_list.append(record.annotations['source'])
               keyw_list.append(record.annotations['keywords'])
               taxo_list.append(record.annotations['taxonomy'])
               desc_list.append(record.description)
               cdsnum+=1
           except:
                   print("Error ",record.description)



#%%
#List of dictionaties to the numpy array
aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
nseqs = len(sequence_list)

percents=sc.zeros((nseqs,20))

for i in range(nseqs):
    percdict = sequence_list[i]
    for an in range(20):
        percents[i,an]= percdict[ aas[an] ]
#plotting percentage heatmap
plt.imshow(percents[:200, :], cmap="hot")


#%%
#PCA
meanX=percents.mean(axis=0)
norm=sc.zeros((percents.shape))
for i in range(percents.shape[0]):
    norm[i,:]=percents[i,:]-meanX
mm = pl.dot(norm.T,norm)
v,e = pl.eig(mm) #eigenvalues, eigenvectors
ev0 = pl.dot(percents,e[:,0])
ev1 = pl.dot(percents,e[:,1])
ev2 = pl.dot(percents,e[:,2])
plt.figure()
plt.plot(ev0,ev1,'.')
plt.figure()
plt.plot(ev0,ev2,'.')
plt.figure()
plt.plot(ev1,ev2,'.')
#plot(percents[:,19],percents[:,18])
#pplot(percents[:,19],percents[:,18],'.')

#%%

plt.figure()
plt.plot(ev0,ev1,'.')
for i in range(len(ev0)):
    plt.text(ev0[i],ev1[i],sources_list[i],fontsize=6)
    

#%%
    
plt.figure()
plt.plot(ev1,ev2,'.')
for i in range(len(ev0)):
    plt.text(ev1[i],ev2[i],desc_list[i],fontsize=6)

#%% 
#Tryptophan
plt.plot(percents[:,18], ".")
#%%

tryptophan=sorted(percents[:,18])
arguments=np.argsort(percents[:,18])
plt.plot(tryptophan, ".")


lowesttryptophan=int(arguments[20:])
highesttryptophan=int(arguments[:20])
#%%
lowesttryptophan_names=[]
for i in range(lowesttryptophan):
    lowesttryptophan_names.append(desc_list[i])

#%%

highesttryptophan_names=[]
for i in range(highesttryptophan):
    highesttryptophan_names.append(desc_list[i])


#%%
def makemaskandplothistogram(desc_list, percen_list, kewd):
    mask=makemasktax(desc_list,kewd)
    mask_perc=[]
    for i in range(len(mask)):
          if (mask[i]):
                 mask_perc.append(percen_list[i,18])
    print(mask_perc)
    plt.hist(mask_perc, bins=int(len(mask_perc)/10))
    plt.hist(mask_perc, bins=int(len(mask_perc)/5))
    plt.hist(mask_perc, bins=int(len(mask_perc)/1))
#%%


makemaskandplothistogram(taxo_list, percents, "aves")

#%%

def makemasktax(alist,keyw):
    mask=sc.zeros(len(alist),bool)
    for i in range(len(alist)):
        arrtostring=','.join(alist[i])
        mask[i]= (arrtostring.upper().find(keyw.upper()) >=0)
    return mask


#%%

def makemaskandplotcdf(desc_list, percen_list, kewd):
    mask=makemask(desc_list,kewd)
    mask_perc=[]
    for i in range(len(mask)):
          if (mask[i]):
                 mask_perc.append(percen_list[i,18])
    sns.kdeplot(mask_perc, shade=True, cumulative=True)
    



#%%
import seaborn as sns
def makemaskandplotcdftax(desc_list, percen_list, kewd):
    mask=makemasktax(desc_list,kewd)
    mask_perc=[]
    for i in range(len(mask)):
          if (mask[i]):
                 mask_perc.append(percen_list[i,18])
    sns.kdeplot(mask_perc, shade=True, cumulative=True)


#%%
makemaskandplotcdftax(taxo_list, percents,"aves")
makemaskandplotcdftax(taxo_list, percents,"rodentia")
makemaskandplotcdftax(taxo_list, percents,"cetacea")

makemaskandplotcdf(desc_list, percents,"alpha")
makemaskandplotcdf(desc_list, percents,"beta")
makemaskandplotcdf(desc_list, percents,"zeta")






plt.title('CDFs of tryptophan percentages in proteins filtered by a keyword')
plt.xlabel('tryptophan percentage')
#%%

