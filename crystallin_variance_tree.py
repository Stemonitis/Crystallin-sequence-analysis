# -*- coding: utf-8 -*-
"""
Spyder Editor
Crystallin variance tree
This is a temporary script file.
"""
#%%

from ete3 import NCBITaxa
from ete3 import Tree

#To prevent console from dying
import os
os.environ['QT_QPA_PLATFORM']='offscreen'


#%%
# Loads a tree structure from a newick string. The returned variable ’t’ is the root node for the tree.
t = Tree("(A:1,(B:1,(E:1,D:1):0.5):0.5);" )
##then to view use
t.render("%%inline")

#%%
 Load a tree structure from a newick file.
t = Tree("genes_tree.nh")
#%%
# You can also specify the newick format. For instance, for named internal nodes we will use format 1.
t = Tree("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", format=1)
#%%
# Loads a tree with internal node names
t = Tree("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", format=1)
#%%
#adddint features (e.g. percentages)

from ete3 import Tree
t =  Tree("((((((4, e), i), o),h), u), ((3, 4), (i, june)));")
# we annotate the tree using external data
colors = {"a":"red", "e":"green", "i":"yellow",
          "o":"black", "u":"purple", "4":"green",
          "3":"yellow", "1":"white", "5":"red",
          "june":"yellow"}
for leaf in t:
    leaf.add_features(color=colors.get(leaf.name, "none"))
print(t.get_ascii(attributes=["name", "color"], show_internal=False))


print("Green-yellow clusters:")
# And obtain clusters exclusively green and yellow
for node in t.get_monophyletic(values=["green", "yellow"], target_attr="color"):
   print(node.get_ascii(attributes=["color", "name"], show_internal=False))
   
   
#%%
   
#finding and saving nodes by their names 
   
C= t&"C"
H= t&"H"
I= t&"I"
   
   
   
#%%
   
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()

#%%

#fetch information
taxid2name = ncbi.get_taxid_translator([9606, 9443])
print(taxid2name)

name2taxid = ncbi.get_name_translator(['Homo sapiens', 'primates'])
print(name2taxid)


# when the same name points to several taxa, all taxids are returned
name2taxid = ncbi.get_name_translator(['Bacteria'])
print(name2taxid)


   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
