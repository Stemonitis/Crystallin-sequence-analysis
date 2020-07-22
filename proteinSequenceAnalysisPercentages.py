# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 14:39:55 2020

@author: mysh
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 22:43:27 2020
Analysis of percentages for protein sequences @author: mysh
"""

#%%
#Importing libraries...
from Bio.SeqIO import parse 
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq 
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
import matplotlib.pyplot as plt
#%%

'''This code is to be utilized in case you already have the list of !!!protein!!! sequences in .gp format and you want to calculate 
percentages of each aminoacid for every sequence''' 


#Specify the path to the file containing sequences
file=open("//home/mysh/colours/scripts/Crystallins/genbankfiles/sequence.gp")
#Creating the list of dictionaries with percentages
sequence_list=[]
for record in  parse(file, "genbank"): #using the parser for .gp genbank file to retrieve and analyze sequences
   current_sequence=str(record.seq) #retrieving the sequence from the record class and transforming it into the string for the subsequent analysis
   analysed_seq = ProteinAnalysis(current_sequence) #creating another fucking class ProteinAnalysis
   sequence_list.append(analysed_seq.get_amino_acids_percent()) #invoking method on this class, it returns a dictionary, we store it in the list
   sequence_list[len(sequence_list)-1]["name"]=str(record.name) #attach the name  of the sequence as on of the keys in the dictionary for every sequence
   sequence_list[-1]["source"]=record.annotations['source']
   #print("Id: %s" % record.id)
   #print("Name: %s" % record.name)
   #print("Description: %s" % record.description)
   #print("Annotations: %s" % record.annotations)
   #print("Sequence Data: %s" % record.seq)
   #print("Sequence Alphabet: %s" % record.seq.alphabet)

#As a result we have the list of dictionaries. Each of this dictionaries corresponds to a sequence.
#Each dictionary has a name of a sequence as one of the key-value pairs and aminoacid percentages 
#as all the the rest 20 of the key-value pairs.

#To count the number of sequences:
number_of_sequences=0
for i in sequence_list:
    number_of_sequences+=1

#To convert the list of dictionaties to the numpy array
sequence_array=np.zeros((number_of_sequences,20))
name_array=[]
acid_array=[]
sequence_number=0
for sequence_dictionary in sequence_list:
    acid_number=0
    for acid in sequence_dictionary:
        if (acid!="name" and acid!="source"):
            sequence_array[sequence_number,acid_number]=sequence_dictionary[acid]
            acid_array.append(acid)
            acid_number+=1
        else:
            name_array.append(sequence_dictionary["name"])
    sequence_number+=1
acid_array=acid_array[:20]
#Now we can further analyze percentages in the numpy array
plt.imshow(sequence_array)
#%%
# To save to npy file
np.save('/home/mysh/colours/crystllines/crystallin_seqs.npy', sequence_array)
np.save('/home/mysh/colours/crystllines/name_array.npy', name_array)
np.save('/home/mysh/colours/crystllines/acid_array.npy', acid_array)