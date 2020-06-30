# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#%%
from ete3 import Tree, TreeStyle, PhyloTree, NCBITaxa, TextFace
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()
#To prevent console from dying while rendering the trees 
import os
os.environ['QT_QPA_PLATFORM']='offscreen'


#%%
taxid_numbers_list=sum(list(ncbi.get_name_translator(taxid_list).values()),[])
#get the tree topology using current ncbi system, make sure you keep all the intermediate nodes
tree = ncbi.get_topology(taxid_numbers_list, intermediate_nodes=True)
#save taxids in the separate attributes and sci_names as actual node names, so they are properly displayed
for node in tree.traverse():
    node.taxid=node.name
    node.name=node.sci_name
    if node.children:
       node.add_face(TextFace(node.sci_name), column=0, position = "branch-right")
    
    
#%%
style = TreeStyle()
#style.mode = "c" # draw tree in circular mode
#style.scale = 20
style.show_leaf_name = False
#uncomment for the 
#style.arc_start = -180 # 0 degrees = 3 o'clock
#style.arc_span = 180
#render the topolgy of species in the tree
tree.render("%%inline", w=800, units="mm", tree_style=style)

#%%
#Adding sequences to the general topology tree

for i in l:
    leaf=tree.search_nodes(taxid=str(i["ncbi_num"][0]))
    leaf.add_child(name=i["name"], ncbi_num=i["ncbi_num"], W=i["W"], N=i["N"], I=i["I"], L=i["L"], D=i["D"], Q=i["Q"],S=i["S"], G=i["G"],V=i["V"], H=i["H"], T=i["T"],A=i["A"], E=i["E"],F=i["F"],K=i["K"],P=i["P"],C=i["C"],Y=i["Y"],R=i["R"],M=i["M"])




