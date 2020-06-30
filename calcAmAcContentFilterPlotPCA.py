
"""
Read, translate, calculate percentages and filter crystallines by keywords in names, sources, taxonomy 
keywords list.
Plot PCA 
"""

#%%
#Importing libraries...
from Bio.SeqIO import parse 
#from Bio.SeqRecord import SeqRecord 
#from Bio.Seq import Seq 
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import scipy.linalg as sc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
#%%
"""What if we want to filter any of the results from the file?"""
#find all elements in a list that match keyw, case INsensitive 
def makemask(alist,keyw):
    mask=np.zeros(len(alist),bool)
    for i in range(len(alist)):
        mask[i]= ( alist[i].upper().find(keyw.upper()) >=0 )
    return mask
    
file=open("/home/mysh/colours/scripts/Crystallins/txt/allgenbank.txt")
file1=open("/home/mysh/colours/scripts/Crystallins/txt/genebankwhalesdolphins.txt")

sequence_list=[]
names_list=[]
sources_list = []
desc_list = []
taxo_list = []
keyw_list = []
for record in  parse(file, "genbank"):
   cdsnum=0
   for feat in record.features:
       if(feat.type == 'CDS'):
           try: 
               prot=feat.translate(record.seq)
               hastrans=True
               analysed_seq = ProteinAnalysis(str(prot)) #creating another fucking class ProteinAnalysis
               sequence_list.append(analysed_seq.get_amino_acids_percent()) #invoking method on this class, it returns a dictionary, we store it in the list
               names_list.append(str(record.name)+ "_CDS#" + str(cdsnum))                    
               sources_list.append(record.annotations['source'])
               keyw_list.append(record.annotations['keywords'])
               taxo_list.append(record.annotations['taxonomy'])
               desc_list.append(record.description)
               cdsnum+=1
           except:
                   print("Error ",record.description) #in case there is a stop codon somewhere there
                   
                   
#List of dictionaties to the numpy array
aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
nseqs = len(sequence_list)

percents=np.zeros((nseqs,20))

for i in range(nseqs):
    percdict = sequence_list[i]
    for an in range(20):
        percents[i,an]= percdict[ aas[an] ]
#%%
plt.semilogy(range(len(v)),np.sort(v), "o", label="with tryptophan")
plt.semilogy(range(len(v1)),np.sort(v1), "o", label="without tryptophan")
plt.legend()
plt.title("Eigenvalues")
#%%
#PCA
x = np.delete(percents, (18), axis=1)
mm = np.dot(x.T,x)
v,e = sc.eig(mm)
ev0 = np.dot(x,e[:,0])
ev1 = np.dot(x,e[:,1])
ev2 = np.dot(x,e[:,2])
#%%
#PCA
mm = np.dot(percents.T,percents)
v1,e = sc.eig(mm)
ev0 = np.dot(x,e[:,0])
ev1 = np.dot(x,e[:,1])
ev2 = np.dot(x,e[:,2])
#%%
# Create 2x2 sub plots
gs = gridspec.GridSpec(2, 2)
pl.title("PCA of protein aminoacid contents in whales&dolphins ")
pl.figure()
ax = pl.subplot(gs[0, 0]) # row 0, col 0
pl.plot(ev0,ev1,'.')
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 2')
ax = pl.subplot(gs[0, 1]) # row 0, col 1
plt.plot(ev0,ev2,'.')
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 3')
ax = pl.subplot(gs[1, 0]) # row 1, span all columns
plt.plot(ev1,ev2,'.')
plt.xlabel('PCA component 2')
plt.ylabel('PCA component 3')

#%%
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(ev0,ev1,ev2, marker='o')
plt.show()
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 2')
#%%
#Plot with names
plt.figure()
plt.plot(ev2,ev1,'.')
for i in range(len(ev0)):
    plt.text(ev0[i],ev1[i],sources_list[i],fontsize=6)

#%%
#Plot with description
plt.figure()
plt.plot(ev1,ev2,'.')
for i in range(len(ev0)):
    plt.text(ev1[i],ev2[i],desc_list[i],fontsize=6)
    
#%%
#Plot colors with types of crystallines
alpha_mask=makemask(desc_list, "alpha")
beta_mask=makemask(desc_list, "beta")
gamma_mask=makemask(desc_list, "gamma")
delta_mask=makemask(desc_list, "delta")
mu_mask=makemask(desc_list, "mu")
lambda_mask=makemask(desc_list, "lambda")
zeta_mask=makemask(desc_list, "zeta")
J1_mask=makemask(desc_list, "J1")

beta_gamma=beta_mask*gamma_mask
complete_mask=makemask(desc_list, "complete")

plt.figure()
plt.plot(ev1,ev2,'.')
plt.plot(ev1*alpha_mask,ev2*alpha_mask, "r.", label="alpha")
plt.plot(ev1*beta_mask,ev2*beta_mask, "g.", color="orchid", label="beta")
plt.plot(ev1*beta_gamma,ev2*beta_gamma, "m.", label="beta and gamma")
plt.plot(ev1*zeta_mask,ev2*zeta_mask, "c.", label="zeta")
plt.plot(ev1*mu_mask,ev2*mu_mask, ".",color='pink', label="mu")
plt.plot(ev1*delta_mask,ev2*delta_mask, ".",color='orange', label="delta")
plt.plot(ev1*lambda_mask,ev2*lambda_mask, ".",color='blueviolet', label="lambda")
plt.plot(ev1*complete_mask,ev2*complete_mask, ".",color='blue', label="complete sequence")
plt.plot(ev1*J1_mask,ev2*J1_mask, ".",color='lightsteelblue', label="J1-")
plt.plot(ev1*gamma_mask,ev2*gamma_mask, ".", color="gold", label="gamma")

plt.legend(shadow=True, ncol=1, fontsize="x-small", markerscale=3)
plt.title("PC 2 and 3 with clustering by the type of the crystallin")
plt.xlabel('PCA component 2')
plt.ylabel('PCA component 3')
#%%
#Plot colors with types of crystallines in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(ev0,ev1,ev2, marker='.', color="blue")
plt.show()
alpha_mask=makemask(desc_list, "alpha")
beta_mask=makemask(desc_list, "beta")
gamma_mask=makemask(desc_list, "gamma")
delta_mask=makemask(desc_list, "delta")
mu_mask=makemask(desc_list, "mu")
lambda_mask=makemask(desc_list, "lambda")
zeta_mask=makemask(desc_list, "zeta")
J1_mask=makemask(desc_list, "J1")
beta_gamma=beta_mask*gamma_mask
complete_mask=makemask(desc_list, "complete")

plt.plot(ev0*alpha_mask,ev1*alpha_mask, ev2*alpha_mask, "r.", label="alpha")
plt.plot(ev0*gamma_mask,ev1*gamma_mask, ev2*gamma_mask, ".", color="gold", label="gamma")
plt.plot(ev0*beta_mask,ev1*beta_mask, ev2*beta_mask, "g.", color="orchid", label="beta")
plt.plot(ev0*beta_gamma,ev1*beta_gamma, ev2*beta_gamma, "m.", label="beta and gamma")
plt.plot(ev0*zeta_mask,ev1*zeta_mask, ev2*zeta_mask, "c.", label="zeta")
plt.plot(ev0*mu_mask,ev1*mu_mask, ev2*mu_mask, ".",color='pink', label="mu")
plt.plot(ev0*delta_mask,ev1*delta_mask,ev2*delta_mask, ".",color='orange', label="delta")
plt.plot(ev0*lambda_mask,ev1*lambda_mask, ev2*lambda_mask, ".",color='blueviolet', label="lambda")
plt.plot(ev0*complete_mask,ev1*complete_mask, ev2*complete_mask, ".",color='blue', label="complete sequence")
plt.plot(ev0*J1_mask,ev1*J1_mask, ev2*J1_mask, ".",color='lightsteelblue', label="J1-")
plt.legend(shadow=True, ncol=1, fontsize="x-small", markerscale=3)
plt.title("PC 1 and 3 with clustering by the type of the crystallin")
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 2')
#%%
#Plot colors with types of crystallines
alpha_mask=makemask(desc_list, "alpha")
beta_mask=makemask(desc_list, "beta")
gamma_mask=makemask(desc_list, "gamma")
delta_mask=makemask(desc_list, "delta")
mu_mask=makemask(desc_list, "mu")
lambda_mask=makemask(desc_list, "lambda")
zeta_mask=makemask(desc_list, "zeta")
J1_mask=makemask(desc_list, "J1")

beta_gamma=beta_mask*gamma_mask
complete_mask=makemask(desc_list, "complete")

plt.figure()
plt.plot(ev0,ev2,'.')
plt.plot(ev0*alpha_mask,ev2*alpha_mask, "r.", label="alpha")
plt.plot(ev0*gamma_mask,ev2*gamma_mask, ".", color="gold", label="gamma")
plt.plot(ev0*beta_mask,ev2*beta_mask, "g.", color="orchid", label="beta")
plt.plot(ev0*beta_gamma,ev2*beta_gamma, "m.", label="beta and gamma")
plt.plot(ev0*zeta_mask,ev2*zeta_mask, "c.", label="zeta")
plt.plot(ev0*mu_mask,ev2*mu_mask, ".",color='pink', label="mu")
plt.plot(ev0*delta_mask,ev2*delta_mask, ".",color='orange', label="delta")
plt.plot(ev0*lambda_mask,ev2*lambda_mask, ".",color='blueviolet', label="lambda")
plt.plot(ev0*complete_mask,ev2*complete_mask, ".",color='blue', label="complete sequence")
plt.plot(ev0*J1_mask,ev2*J1_mask, ".",color='lightsteelblue', label="J1-")
plt.legend(shadow=True, ncol=1, fontsize="x-small", markerscale=3)
plt.title("PC 1 and 3 with clustering by the type of the crystallin")
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 3')
#%%
#Plot colors with taxonomic groups
#case senstive
def makemasktax(taxo_list,keyw):
    mask=np.zeros(len(taxo_list),bool)
    for i in range(len(taxo_list)):
        if (any(l==keyw for l in taxo_list[i])):
            mask[i]= True
    return mask

aves_mask=makemasktax(taxo_list, "Aves")
dol_mask=makemasktax(taxo_list, "Tursiops")
mam_mask=makemasktax(taxo_list, "Mammalia")
rod_mask=makemasktax(taxo_list, "Rodentia")
fish_mask=makemasktax(taxo_list, "Actinopterygii")
rep_mask=makemasktax(taxo_list, "Reptilia")
frog_mask=makemasktax(taxo_list, "Amphibia")
nem_mask=makemasktax(taxo_list, "Nematoda")
ar_mask=makemasktax(taxo_list, "Archelosauria")
plant_mask=makemasktax(taxo_list, "Viridiplantae")
cet_mask=makemasktax(taxo_list, "Cetacea")



plt.figure()
plt.plot(ev0,ev1,'.')
plt.plot(ev0*mam_mask,ev1*mam_mask, ".",color="orchid", label="Mammalia")
plt.plot(ev0*ar_mask,ev1*ar_mask, ".",color="lightsteelblue", label="Archelosauria")
plt.plot(ev0*aves_mask,ev1*aves_mask, ".",color="gold", label="Aves")
plt.plot(ev0*rod_mask,ev1*rod_mask, "m.", label="Rodentia")
plt.plot(ev0*fish_mask,ev1*fish_mask, "c.", label="Actinopterygii")
plt.plot(ev0*rep_mask,ev1*rep_mask, ".",color="pink", label="Reptilia")
plt.plot(ev0*frog_mask,ev1*frog_mask, ".",color="orange", label="Amphibia")
plt.plot(ev0*nem_mask,ev1*nem_mask, ".",color="blueviolet", label="Nematoda")
plt.plot(ev0*plant_mask,ev1*plant_mask, ".",color="darkgreen", label="Viridaplantae")
plt.plot(ev0*cet_mask,ev1*cet_mask, ".",color='lightsteelblue', label="Cetacea")
plt.plot(ev0*dol_mask,ev1*dol_mask, "r.", label="Tursiops")



plt.legend(shadow=True, ncol=1, fontsize="x-small", markerscale=3)
plt.title("PC 1 and 2 with clustering by taxonomic group")
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 1')


    
#%%
    
#Plot colors with ttaxonomy in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(ev0,ev1,ev2, marker='.')
plt.show()
aves_mask=makemasktax(taxo_list, "Aves")
dol_mask=makemasktax(taxo_list, "Tursiops")
mam_mask=makemasktax(taxo_list, "Mammalia")
rod_mask=makemasktax(taxo_list, "Rodentia")
fish_mask=makemasktax(taxo_list, "Actinopterygii")
rep_mask=makemasktax(taxo_list, "Reptilia")
frog_mask=makemasktax(taxo_list, "Amphibia")
nem_mask=makemasktax(taxo_list, "Nematoda")
ar_mask=makemasktax(taxo_list, "Archelosauria")
plant_mask=makemasktax(taxo_list, "Viridiplantae")
cet_mask=makemasktax(taxo_list, "Cetacea")

plt.plot(ev0*aves_mask,ev1*aves_mask, ev2*aves_mask, "r.", label="Aves")
plt.plot(ev0*dol_mask,ev1*dol_mask, ev2*dol_mask, ".", color="gold", label="Tursiops")
plt.plot(ev0*mam_mask,ev1*mam_mask, ev2*mam_mask, "g.", color="orchid", label="Mammalia")
plt.plot(ev0*rod_mask,ev1*rod_mask, ev2*rod_mask, "m.", label="Rodentia")
plt.plot(ev0*fish_mask,ev1*fish_mask, ev2*fish_mask, "c.", label="Actinopterygii")
plt.plot(ev0*rep_mask,ev1*rep_mask, ev2*rep_mask, ".",color='pink', label="Reptilia")
plt.plot(ev0*frog_mask,ev1*frog_mask,ev2*frog_mask, ".",color='orange', label="Amphibia")
plt.plot(ev0*nem_mask,ev1*nem_mask, ev2*nem_mask, ".",color='blueviolet', label="Nematoda")
plt.plot(ev0*ar_mask,ev1*ar_mask, ev2*ar_mask, ".",color='blue', label="Archelosauria")
plt.plot(ev0*plant_mask,ev1*plant_mask, ev2*plant_mask, ".",color='lightsteelblue', label="Viridiplantae")
plt.plot(ev0*cet_mask,ev1*cet_mask, ev2*cet_mask, ".",color='lightsteelblue', label="Cetacea")

plt.legend(shadow=True, ncol=1, fontsize="x-small", markerscale=3)
plt.title("PC 1 and 2 with clustering by the taxonomy group")
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 2')

#%%
#Tryptophan content
tr_perc=percents[:,18]
c=tr_perc


# Create 2x2 sub plots
gs = gridspec.GridSpec(2, 2)
plt.title("PC 1 and 2 with clustering by the tryptophan content")

pl.figure()
ax = pl.subplot(gs[0, 0]) # row 0, col 0
plt.scatter(ev0,ev1,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 2')
ax = pl.subplot(gs[0, 1]) # row 0, col 1
plt.scatter(ev0,ev2,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 3')
ax = pl.subplot(gs[1, 0]) # row 1, span all columns
plt.scatter(ev1,ev2,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.xlabel('PCA component 2')
plt.ylabel('PCA component 3')
plt.colorbar()

#%%
#Tryptophan content 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(ev0,ev1,ev2, c=tr_perc,cmap = 'seismic',s=50, alpha=0.5 )
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 2')
plt.title("PCA tryptophan percentage")

plt.colorbar()
plt.show()

#%%
all_mask=[not i for i in (aves_mask+cet_mask+mam_mask+rod_mask+fish_mask+rep_mask+frog_mask+nem_mask+ar_mask)]
#%%
for i in range(len(ev0)):
    if all_mask[i]:
        plt.plot(ev0[i], ev2[i])
        plt.text(ev0[i],ev2[i],taxo_list[i],fontsize=6)
        #%%
#Plot colors with taxonomic groups
#case senstive
def makemasktax(taxo_list,keyw):
    mask=np.zeros(len(taxo_list),bool)
    for i in range(len(taxo_list)):
        if (any(l==keyw for l in taxo_list[i])):
            mask[i]= True
    return mask

aves_mask=makemasktax(taxo_list, "Sousa")
dol_mask=makemasktax(taxo_list, "Tursiops")
rod_mask=makemasktax(taxo_list, "Monodon")
phys_mask=makemasktax(taxo_list, "Physeter")
del_mask=makemasktax(taxo_list, "Delphinapterus")
bal_mask=makemasktax(taxo_list, "Balaenoptera")



plt.figure()
plt.plot(ev1,ev2,'.')
plt.plot(ev1*aves_mask,ev2*aves_mask, ".",color="gold", label="Sousa")
plt.plot(ev1*rod_mask,ev2*rod_mask, "m.", label="Monodon")
plt.plot(ev1*phys_mask,ev2*phys_mask, "c.", label="Physeter")
plt.plot(ev1*dol_mask,ev2*dol_mask, "r.", label="Tursiops")
plt.plot(ev1*del_mask,ev2*del_mask,".", color="lightsteelblue", label="Delphinapterus")
plt.plot(ev1*bal_mask,ev2*bal_mask,".g", label="Balaenoptera")



plt.legend(shadow=True, ncol=1, fontsize="x-small", markerscale=3)
plt.title("PC 2 and 3 with clustering by taxonomic group")
plt.xlabel('PCA component 2')
plt.ylabel('PCA component 3')

#%%
#Plot colors with ttaxonomy in 3D for whales and dolphins
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(ev0,ev1,ev2, marker='.')
plt.show()
aves_mask=makemasktax(taxo_list, "Sousa")
dol_mask=makemasktax(taxo_list, "Tursiops")
rod_mask=makemasktax(taxo_list, "Monodon")
phys_mask=makemasktax(taxo_list, "Physeter")
del_mask=makemasktax(taxo_list, "Delphinapterus")
bal_mask=makemasktax(taxo_list, "Balaenoptera")

plt.plot(ev0*aves_mask,ev1*aves_mask, ev2*aves_mask, "r.", label="Sousa")
plt.plot(ev0*dol_mask,ev1*dol_mask, ev2*dol_mask, ".", color="black", label="Tursiops")
plt.plot(ev0*rod_mask,ev1*rod_mask, ev2*rod_mask, "m.", label="Monodon")
plt.plot(ev0*phys_mask,ev1*phys_mask, ev2*phys_mask, ".",color='pink', label="Physeter")
plt.plot(ev0*del_mask,ev1*del_mask,ev2*del_mask, ".",color='orange', label="Delphinapterus")
plt.plot(ev0*bal_mask,ev1*bal_mask, ev2*bal_mask, ".",color='blueviolet', label="Balaenoptera")

plt.legend(shadow=True, ncol=1, fontsize="x-small", markerscale=3)
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 2')
#%%

#Alpha-crystallines taxonomy group tryptophan content 3D

aves_mask=makemasktax(taxo_list, "Aves")
dol_mask=makemasktax(taxo_list, "Tursiops")
mam_mask=makemasktax(taxo_list, "Mammalia")
rod_mask=makemasktax(taxo_list, "Rodentia")
fish_mask=makemasktax(taxo_list, "Actinopterygii")
rep_mask=makemasktax(taxo_list, "Reptilia")
frog_mask=makemasktax(taxo_list, "Amphibia")
nem_mask=makemasktax(taxo_list, "Nematoda")
ar_mask=makemasktax(taxo_list, "Archelosauria")
plant_mask=makemasktax(taxo_list, "Viridiplantae")
cet_mask=makemasktax(taxo_list, "Cetacea")
alpha_mask=makemask(desc_list, "alpha")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.show()

for i in range(len(ev0)):
    if alpha_mask[i]:
        ax.scatter(ev0[i],ev1[i],ev2[i], marker='.')
        ax.scatter(ev0[i]*aves_mask[i],ev1[i]*aves_mask[i], ev2[i]*aves_mask[i], "r.", label="Aves")
        ax.scatter(ev0[i]*dol_mask[i],ev1[i]*dol_mask[i], ev2[i]*dol_mask[i], ".", color="gold", label="Tursiops")
        ax.scatter(ev0[i]*mam_mask[i],ev1[i]*mam_mask[i], ev2[i]*mam_mask[i], "g.", color="orchid", label="Mammalia")
        ax.scatter(ev0[i]*rod_mask[i],ev1[i]*rod_mask[i], ev2[i]*rod_mask[i], "m.", label="Rodentia")
        ax.scatter(ev0[i]*fish_mask[i],ev1[i]*fish_mask[i], ev2[i]*fish_mask[i], "c.", label="Actinopterygii")
        ax.scatter(ev0[i]*rep_mask[i],ev1[i]*rep_mask[i], ev2[i]*rep_mask[i], ".",color='pink', label="Reptilia")
        ax.scatter(ev0[i]*frog_mask[i],ev1[i]*frog_mask[i],ev2[i]*frog_mask[i], ".",color='orange', label="Amphibia")
        ax.scatter(ev0[i]*nem_mask[i],ev1[i]*nem_mask[i], ev2[i]*nem_mask[i], ".",color='blueviolet', label="Nematoda")
        ax.scatter(ev0[i]*ar_mask[i],ev1[i]*ar_mask[i], ev2[i]*ar_mask[i], ".",color='blue', label="Archelosauria")
        ax.scatter(ev0[i]*plant_mask[i],ev1[i]*plant_mask[i], ev2[i]*plant_mask[i], ".",color='lightsteelblue', label="Viridiplantae")
        ax.scatter(ev0[i]*cet_mask[i],ev1[i]*cet_mask[i], ev2[i]*cet_mask[i], ".",color='lightsteelblue', label="Cetacea")
        
#plt.legend(shadow=True, ncol=1, fontsize="x-small", markerscale=3)
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 2')
     
#%%
#Alpha-crystallines tryptophan content
#Tryptophan content
tr_perc=percents[:,18]
c=tr_perc
ev0=ev1
ev1=ev2
alpha_mask=makemask(desc_list, "alpha")
beta_mask=makemask(desc_list, "beta")
gamma_mask=makemask(desc_list, "gamma")
delta_mask=makemask(desc_list, "delta")
mu_mask=makemask(desc_list, "mu")
lambda_mask=makemask(desc_list, "lambda")
zeta_mask=makemask(desc_list, "zeta")
J1_mask=makemask(desc_list, "J1")
beta_gamma=beta_mask*gamma_mask
complete_mask=makemask(desc_list, "complete")
ax = pl.subplot(gs[0, 0]) # row 0, col 0
plt.scatter(ev0*alpha_mask,ev1*alpha_mask,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.title("Alpha")
ax = pl.subplot(gs[0, 1]) # row 0, col 0
plt.scatter(ev0*beta_mask,ev1*beta_mask,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.title("Beta")
ax = pl.subplot(gs[0, 2]) # row 0, col 0
plt.scatter(ev0*gamma_mask,ev1*gamma_mask,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.title("Gamma")
ax = pl.subplot(gs[1, 0]) # row 0, col 0
plt.scatter(ev0*delta_mask,ev1*delta_mask,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.title("Delta")
ax = pl.subplot(gs[1, 1]) # row 0, col 0
plt.scatter(ev0*mu_mask,ev1*mu_mask,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.title("Mu")
ax = pl.subplot(gs[1, 2]) # row 0, col 0
plt.scatter(ev0*lambda_mask,ev1*lambda_mask,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.title("Lambda")
ax = pl.subplot(gs[2, 0]) # row 0, col 0
plt.scatter(ev0*zeta_mask,ev1*zeta_mask,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.title("Zeta")
ax = pl.subplot(gs[2, 1]) # row 0, col 0
plt.scatter(ev0*beta_gamma,ev1*beta_gamma,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.title("Beta and gamma")
ax = pl.subplot(gs[2, 1]) # row 0, col 0
plt.scatter(ev0*J1_mask,ev1*J1_mask,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.title("J1")

plt.colorbar()


#%%




# Create 2x2 sub plots
gs = gridspec.GridSpec(2, 2)
plt.title("PC 1 and 2 with clustering by the tryptophan content")

pl.figure()
ax = pl.subplot(gs[0, 0]) # row 0, col 0
plt.scatter(ev0,ev1,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 2')
ax = pl.subplot(gs[0, 1]) # row 0, col 1
plt.scatter(ev0,ev2,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 3')
ax = pl.subplot(gs[1, 0]) # row 1, span all columns
plt.scatter(ev1,ev2,c=tr_perc,cmap = 'seismic',s=50, alpha=0.5)
plt.xlabel('PCA component 2')
plt.ylabel('PCA component 3')
plt.colorbar()











