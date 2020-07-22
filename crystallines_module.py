# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 08:20:25 2020

@author: mysh
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
import statistics

<<<<<<< HEAD
from ete3 import Tree, TreeStyle, PhyloTree, NCBITaxa, TextFace, NodeStyle
=======
from ete3 import Tree, TreeStyle, PhyloTree, NCBITaxa, TextFace
>>>>>>> a1e6c83dc21f1ab374ad7d39a11391f872d50b87
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()
#To prevent console from dying while rendering the trees 
import os
os.environ['QT_QPA_PLATFORM']='offscreen'
#%%

"""What if we want to filter any of the results from the file?"""
#find all elements in a list that match keyw, case INsensitive 
def makemask(alist,keyw):
    mask=np.zeros(len(alist),bool)
    for i in range(len(alist)):
        mask[i]= ( alist[i].upper().find(keyw.upper()) >=0 )
    return mask
    
def makemasknestedlist(taxo_list,keyw):
    mask=np.zeros(len(taxo_list),bool)
    for i in range(len(taxo_list)):
        if (any(l==keyw for l in taxo_list[i])):
            mask[i]= True
    return mask

def percentages_from_proteins(path):
    file=open(path)
    names_list=[]
    sequence_list=[]
    sources_list = []
    desc_list = []
    taxo_list = []
    keyw_list = []
    taxid_list = []
    for record in  parse(file, "genbank"):
      cdsnum=0
      for feat in record.features:
               prot=record.seq
               analysed_seq = ProteinAnalysis(str(prot)) #creating another class ProteinAnalysis
               sequence_list.append(analysed_seq.get_amino_acids_percent()) #invoking method on this class, it returns a dictionary, we store it in the list
               names_list.append(str(record.name)+ "_CDS#" + str(cdsnum))                    
               sources_list.append(record.annotations['source'])
               keyw_list.append(record.annotations['keywords'])
               taxo_list.append(record.annotations['taxonomy'])
               desc_list.append(record.description)
               taxid_list.append(record.annotations["organism"])
               cdsnum+=1
    #List of dictionaties to the numpy array
    aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    nseqs = len(sequence_list)
    percents=np.zeros((nseqs,20))
    for i in range(nseqs):
        percdict = sequence_list[i]
        for an in range(20):
             percents[i,an]= percdict[ aas[an] ]
    return percents, names_list, sources_list, desc_list, taxo_list, keyw_list, taxid_list, sequence_list

def percentages_from_nucleotides(path):
   file=open(path)
   names_list=[]
   sequence_list=[]
   sources_list = []
   desc_list = []
   taxo_list = []
   keyw_list = []
   taxid_list = []
   all_attr_dict={}

   for record in  parse(file, "genbank"):
    cdsnum=0
    for feat in record.features:
       if(feat.type == 'CDS'):
           try: 
               prot=feat.translate(record.seq)
               analysed_seq = ProteinAnalysis(str(prot)) #creating another class ProteinAnalysis
               sequence_list.append(analysed_seq.get_amino_acids_percent()) #invoking method on this class, it returns a dictionary, we store it in the list
               names_list.append(str(record.name)+ "_CDS#" + str(cdsnum))                    
               sources_list.append(record.annotations['source'])
               keyw_list.append(record.annotations['keywords'])
               taxo_list.append(record.annotations['taxonomy'])
               desc_list.append(record.description)
               taxid_list.append(record.source["organism"])
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
    return percents, names_list, sources_list, desc_list, taxo_list, keyw_list, taxid_list, sequence_list
    
def PCA(percents):
    mm = np.dot(percents.T,percents)
    v1,e = sc.eig(mm)
    ev0 = np.dot(percents,e[:,0])
    ev1 = np.dot(percents,e[:,1])
    ev2 = np.dot(percents,e[:,2])
    return ev0, ev1, ev2, v1
<<<<<<< HEAD
#where n is a list of principal components, don`t start with 0!
def plotPC(percents, n):
    mm = np.dot(percents.T,percents)
    v1,e = sc.eig(mm)
    for i in range(len(n)):
        plt.plot(e[i], ".-", label="Principal components "+str(i+1))
    aa = ['A(Ala)', 'C(Cys)', 'D(Asp)', 'E(Glu)', 'F(Phe)', 'G(Gly)', 'H(His)', 'I(Ile)', 'K(Lys)', 'L(Leu)', 'M(Met)', 'N(Asn)', 'P(Pro)', 'Q(Gln)', 'R(Arg)', 'S(Ser)', 'T(Thr)', 'V(Val)', 'W(Trp)', 'Y(Tyr)']
    ai=list(range(len(aa)))
    plt.xlabel('Amino acid type')
    plt.xticks(ai, aa)
    plt.legend()
    plt.title('Principal components')
=======
>>>>>>> a1e6c83dc21f1ab374ad7d39a11391f872d50b87
#See how absence of one of the aminoacids influences the PCA aa--number of the aa in the list, e.g. tryptophan=18
def PCA_without(percents, aa):
    x = np.delete(percents, (aa), axis=1)
    mm = np.dot(x.T,x)
    v,e = sc.eig(mm)
    ev0 = np.dot(x,e[:,0])
    ev1 = np.dot(x,e[:,1])
    ev2 = np.dot(x,e[:,2])
    return ev0,ev1, ev2,v
    
def plot_eigenvalues(v):
    plt.semilogy(range(len(v)),np.sort(v), "o")
    plt.title("Eigenvalues")

def plot_with_without_tryptophan(percents):
    ev0, ev1, ev2, v=PCA(percents)
    ev0, ev1, ev2, v1=PCA_without(percents, 18)
    plt.semilogy(range(len(v)),np.sort(v), "o", label="with tryptophan")
    plt.semilogy(range(len(v1)),np.sort(v1), "o", label="without tryptophan")
    plt.legend()
    plt.title("Eigenvalues")
    
def plot_PCA_components(percents,title):
    ev0, ev1, ev2, v=PCA(percents)
    gs = gridspec.GridSpec(2, 2)
    pl.title("gs = gridspec.GridSpec(2, 2)")
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
    
def plot_PCA_components3D(percents,title):
    ev0, ev1, ev2, v=PCA(percents)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(ev0,ev1,ev2, marker='o')
    plt.title(title)
    plt.show()
    plt.xlabel('PCA component 1')
    plt.ylabel('PCA component 2')
    
def plot_with_text(percents, txt_list):
    ev0, ev1, ev2, v=PCA(percents)
    plt.figure()
    plt.plot(ev0,ev1,'.')
    for i in range(len(ev0)):
        plt.text(ev0[i],ev1[i],txt_list[i],fontsize=6)
#PC_number tells is a list of two principal components that you want to plot. Possible entries:
#[0,1] -- first and second, [0,2] -- first and third, [1,2] -- second and third.
#cat_group is the string that characterizes the groups name (e.g. "taxonomic group", "type of crystallin", etc.)
#If there are overlapping groups than the latest one will overlap all the previuos ones
def plot_with_mask(percents, categories_list, search_group, search_list, PC_number):
    ev0, ev1, ev2, v=PCA(percents)
    PCs=[ev0,ev1,ev2]
    plt.figure()
    plt.plot(PCs[PC_number[0]],PCs[PC_number[1]],'.')
    for cat in categories_list:
        #check if the search_list is nested
        if type(search_list[0])==str:
            mask=makemask(search_list, cat)
        else:
            mask=makemasknestedlist(search_list,cat)
        plt.plot(PCs[PC_number[0]]*mask,PCs[PC_number[1]]*mask, ".", label=cat)
        plt.legend(shadow=True, ncol=1, fontsize="x-small", markerscale=3)
    plt.title("PC 1 and 2 with clustering by "+search_group)
    plt.xlabel('PCA component '+str(PC_number[0]))
    plt.ylabel('PCA component '+str(PC_number[1]))

def plot_with_mask3D(percents, categories_list, search_group, search_list):
    ev0, ev1, ev2, v=PCA(percents)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(ev0,ev1,ev2, marker='.')
    plt.show()
    for cat in categories_list:
        #check if the search_list is nested
        if type(search_list[0])==str:
            mask=makemask(search_list, cat)
        else:
            mask=makemasknestedlist(search_list,cat)
        plt.plot(ev0*mask,ev1*mask, ev2*mask, ".", label=cat)
        plt.legend(shadow=True, ncol=1, fontsize="x-small", markerscale=3)
    plt.title("PC 1 and 2 with clustering by "+search_group)
    plt.xlabel('PCA component 1')
    plt.ylabel('PCA component 2')
    
<<<<<<< HEAD
def plot_aa_content(percents,aa):
    ev0, ev1, ev2, v=PCA(percents)
    aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    aa_perc=percents[:,aas.index(aa)]
    # Create 2x2 sub plots
    gs = gridspec.GridSpec(2, 2)
    ax = pl.subplot(gs[0, 0]) # row 0, col 0
    plt.scatter(ev0,ev1,c=aa_perc,cmap = 'seismic',s=50, alpha=0.5)
    plt.xlabel('PCA component 1')
    plt.ylabel('PCA component 2')
    ax = pl.subplot(gs[0, 1]) # row 0, col 1
    plt.scatter(ev0,ev2,c=aa_perc,cmap = 'seismic',s=50, alpha=0.5)
    plt.xlabel('PCA component 1')
    plt.ylabel('PCA component 3')
    ax = pl.subplot(gs[1, 0]) # row 1, span all columns
    plt.scatter(ev1,ev2,c=aa_perc,cmap = 'seismic',s=50, alpha=0.5)
=======
def plot_tryptophan_content(percents):
    ev0, ev1, ev2, v=PCA(percents)
    tr_perc=percents[:,18]
    c=tr_perc
    # Create 2x2 sub plots
    gs = gridspec.GridSpec(2, 2)
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
>>>>>>> a1e6c83dc21f1ab374ad7d39a11391f872d50b87
    plt.suptitle("PC 1 and 2 with clustering by the tryptophan content")
    plt.xlabel('PCA component 2')
    plt.ylabel('PCA component 3')
    plt.colorbar()
    
<<<<<<< HEAD
def plot_aa_content3D(percents, aa):
    ev0, ev1, ev2, v=PCA(percents)
    aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    aa_perc=percents[:,aas.index(aa)]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(ev0,ev1,ev2, c=aa_perc,cmap = 'seismic',s=50, alpha=0.5 )
=======
def plot_tryptophan_content3D(percents):
    tr_perc=percents[:,18]
    c=tr_perc
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(ev0,ev1,ev2, c=tr_perc,cmap = 'seismic',s=50, alpha=0.5 )
>>>>>>> a1e6c83dc21f1ab374ad7d39a11391f872d50b87
    plt.xlabel('PCA component 1')
    plt.ylabel('PCA component 2')
    plt.title("PCA tryptophan percentage")
    plt.colorbar()
    plt.show()
def create_dictionary_list(percents, names_list, taxid_list, aa_dic):
    dictionary_list=[]
    for i in range(len(percents)):
        dictionary_list.append({"name":names_list[i], "ncbi_num":sum(list(ncbi.get_name_translator([taxid_list[i]]).values()),[])})
        dictionary_list[i].update(aa_dic[i])
    return dictionary_list

def get_sequence_tree(percents, names_list, taxid_list, aa_dic):
    dictionary_list=create_dictionary_list(percents, names_list, taxid_list, aa_dic)
    taxid_numbers_list=sum(list(ncbi.get_name_translator(taxid_list).values()),[])
    #get the tree topology using current ncbi system, make sure you keep all the intermediate nodes
    tree = ncbi.get_topology(taxid_numbers_list, intermediate_nodes=True)
    #save taxids in the separate attributes and sci_names as actual node names, so they are properly displayed
    for node in tree.traverse():
         node.taxid=node.name
         node.name=node.sci_name
         if node.children:
            node.add_face(TextFace(node.sci_name), column=0, position = "branch-right")
    for i in dictionary_list:
       leaf=tree.search_nodes(taxid=str(i["ncbi_num"][0]))[0]
       leaf.add_child(name=i["name"]).add_features(ncbi_num=i["ncbi_num"], W=i["W"], N=i["N"], I=i["I"], L=i["L"], D=i["D"], Q=i["Q"],S=i["S"], G=i["G"],V=i["V"], H=i["H"], T=i["T"],A=i["A"], E=i["E"],F=i["F"],K=i["K"],P=i["P"],C=i["C"],Y=i["Y"],R=i["R"],M=i["M"])
    return tree

<<<<<<< HEAD
def plot_tree(tree, path):
=======

def plot_tree_full(tree, path):
>>>>>>> a1e6c83dc21f1ab374ad7d39a11391f872d50b87
    style = TreeStyle()
    #style.mode = "c" # draw tree in circular mode
    #style.scale = 50
    style.show_leaf_name = False
    #uncomment for the 
    #style.arc_start = -180 # 0 degrees = 3 o'clock
    #style.arc_span = 180
    #render the topolgy of species in the tree
<<<<<<< HEAD
    tree.render(path, w=700, units="mm", tree_style=style)
  
def rgb(minimum, maximum, value):
    minimum, maximum = float(minimum), float(maximum)
    ratio = 2 * (value-minimum) / (maximum - minimum)
    b = int(max(0, 255*(1 - ratio)))
    r = int(max(0, 255*(ratio - 1)))
    g = 255 - b - r
    return r, g, b
    
def red_gradient(minimum, maximum, value):
    minimum, maximum = float(minimum), float(maximum)
    ratio = 2 * (value-minimum) / (maximum - minimum)
    b = int(max(0, 255*(1 - ratio)))
    r = int(max(0, 255*(ratio - 1)))
    g = 255 - b - r
    return r, 255, 255
    
def rgb_to_hex(rgb):
    return '%02x%02x%02x' % rgb
#%%
def plot_variance(tree,path,aa):
    #cycle through the nodes and add a variance attribute to each node postorder, so we start with the
#children that have W attribute and recursively create this attribute as we are traversing through the tree
    all_means=[]
    for node in tree.traverse(strategy="postorder"):
        if not(node.is_leaf()):
            children=node.children
            variance=[]
            for i in range(len(children)):
                aa_feat = getattr(children[i], aa)
                variance.append(aa_feat)
                # variance requires at least two data points, in case there is only one sequence or node
                if len(variance)==1:
                    node.add_feature(aa+"_var",0)
                    node.add_feature(aa+"", np.mean(variance))
                else:
                    node.add_feature(aa+"_var",statistics.variance(variance))
                    node.add_feature(aa+"", np.mean(variance))
            all_means.append(getattr(node, aa+""))
    #plot without sequences (create new tree with cut leave and )
    new_tree=tree.copy()
    maximum=max(all_means)
    minimum=min(all_means)
    #for node in new_tree.traverse(strategy="preorder"):
    ts = TreeStyle()
    ts.mode="c"
    for node in new_tree.traverse():
        if node.is_leaf():
            node.detach()
        else:
            nstyle = NodeStyle()
            nstyle["bgcolor"]="#"+rgb_to_hex(red_gradient(minimum,maximum,(getattr(node, aa+""))))
            node.add_face(TextFace(getattr(node,aa+"")), column=0, position = "branch-right")
            node.set_style(nstyle)
    new_tree.render(path, w=1400, units="mm", tree_style=ts)
#%%
path="/home/mysh/colours/scripts/Crystallins/genbankfiles/crybb1_mammals.gp"
#path2="/home/mysh/colours/scripts/Crystallins/txt/allgenbank.txt"
percents, names_list, sources_list, desc_list, taxo_list, keyw_list, taxid_list, aa_dic =percentages_from_proteins(path)
#percents, names_list, sources_list, desc_list, taxo_list, keyw_list=percentages_from_proteins(path2)
#%%
tree=get_sequence_tree(percents, names_list, taxid_list, aa_dic)
path="/home/mysh/colours/scripts/Crystallins/tree.png"
#render the tree
plot_tree(tree, path)
=======
    tree.render(path, w=400, units="mm", tree_style=style)
#%% 
    
def plot_variance(tree,path,aa):
    #cycle through the nodes and add a variance attribute to each node
    for node in tree.traverse(strategy="postorder"):
        if !(node.is_leaf):
            children=node.children
            variance=[]
            for i in children:
                variance.append(i.aa)
            node.add_feature[aa+""=statistics.variance(variance)]
    #plot without sequences (create new tree with cut leave and )
    new_tree=tree.copy()
    for node in tree:
        if tree.node_is_leaf():
            
            
            

#%%
path="/home/mysh/colours/scripts/Crystallins/genbankfiles/cryaa_mammals.gp"
#path2="/home/mysh/colours/scripts/Crystallins/txt/allgenbank.txt"
percents, names_list, sources_list, desc_list, taxo_list, keyw_list, taxid_list, aa_dic =percentages_from_proteins(path)
#percents, names_list, sources_list, desc_list, taxo_list, keyw_list=percentages_from_proteins(path2)
tree=get_sequence_tree(percents, names_list, taxid_list, aa_dic)
#%%
path="/home/mysh/colours/scripts/Crystallins/mytree.png"
#render the tree
plot_tree(tree, path)
#%%

#%%
>>>>>>> a1e6c83dc21f1ab374ad7d39a11391f872d50b87
ev0, ev1, ev2, v=PCA(percents)
#ev0, ev1, ev2, v=PCA_without(percents, 18)
#plot_eigenvalues(v)
#plot_with_without_tryptophan(percents)
#plot_PCA_components(percents, "PCA of protein aminoacid contents in")
plot_PCA_components3D(percents, "PCA of protein aminoacid contents in")
plot_with_text(percents,desc_list)
plot_with_mask(percents, ["gamma", "beta"], "type_of_crystallin", desc_list, [0,1])
plot_with_mask(percents, ["Mammalia", "Cetacea"], "taxonomic_group", taxo_list, [0,1])
plot_with_mask3D(percents, ["Mammalia","Cetacea"], "taxonomic_group", taxo_list)
plot_tryptophan_content(percents)
<<<<<<< HEAD
plot_aa_content3D(percents, "W")
#%%


path="/home/mysh/colours/scripts/Crystallins/var_tree.png"
plot_variance(tree,path,"Y")
=======
plot_tryptophan_content3D(percents)

#%%
>>>>>>> a1e6c83dc21f1ab374ad7d39a11391f872d50b87
