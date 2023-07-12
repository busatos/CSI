from Bio import SeqIO, AlignIO
from Bio.SubsMat import MatrixInfo
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.Applications import PhymlCommandline

import dendropy
from dendropy.calculate import treemeasure
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt 
import numpy as np

from random import sample
import sys
import os
import glob
import math
import re 

###a couple things that I need to refer to####
##MW for AAs
my_aa = {
"A": [71.0788, 0],
"R": [156.1875, 0],
"N": [114.1038, 0],
"D": [115.0886, 0],
"C": [103.1388, 0],
"E": [129.1155, 0],
"Q": [128.1307, 0],
"G": [57.0519, 0],
"H": [137.1411, 0],
"I": [113.1594, 0],
"L": [113.1594, 0],
"K": [128.1741, 0],
"M": [131.1926, 0],
"F": [147.1766, 0],
"P": [97.1167, 0],
"S": [87.0782, 0],
"T": [101.1051, 0],
"W": [186.2132, 0],
"Y": [163.176, 0],
"V": [99.1326, 0],
}
##pkvalues
my_pkvals = {
"Nterm": 9.094,
"Cterm": 2.869,
"C": 7.555,
"D": 3.872,
"E": 4.412,
"H": 5.637,
"K": 9.052,
"R": 11.84,
"Y": 10.85,
}

#calculate MW
def mw(seq):
    i = 0
    total_MW = 0
    #my_weights = my_aa
    for i in range(0, len(seq)):
        new_MW = my_aa[str(seq[i])][0]
        total_MW = round(total_MW + new_MW,4)
    return total_MW

#calculate charge based on pH
def charge(seq, pH):
    i = 0
    total_charge = 0
    for i in range (0, len(seq)):
        if seq[i] == "K" or seq[i] == "R" or seq[i] == "H":
            partial_charge = (1/(1 + pow(10,(pH-my_pkvals[str(seq[i])]))))
            total_charge = total_charge + partial_charge
        elif seq[i] == "D" or seq[i] == "E" or seq[i] == "C" or seq[i] == "Y":
            partial_charge = (-1/(1 + pow(10,(my_pkvals[str(seq[i])]-pH))))
            total_charge = total_charge + partial_charge
    n_term_charge = (1/(1 + pow(10,(pH-my_pkvals["Nterm"]))))
    c_term_charge = (-1/(1 + pow(10,(my_pkvals["Cterm"]-pH))))
    total_charge = n_term_charge + total_charge + c_term_charge
    return total_charge

#calculate isoelectric point
def pI(seq):
    test_pH = 7
    diff = 4
    test_charge = charge(seq, test_pH)
    while round(test_charge, 4) != 0:
        if test_charge > 0:
            test_pH = test_pH + diff
            test_charge = charge(seq, test_pH)
            diff = float(diff/2)
        else:
            test_pH = test_pH - diff
            test_charge = charge(seq, test_pH)
            diff = float(diff/2)
    pI_val = round(test_pH, 4)
    return pI_val

#output MW, pI and charge at 7 as a list
def info(seq):
    molw = mw(seq)
    pot = pI(seq)
    ch = round(charge(seq, 7),4)
    results = [molw, pot, ch]
    return results
    
    
##functions to sort anything based on re match of numbers in it

def return_number(text):
    return int(text) if text.isdigit() else text
        
def natural_keys(text):
    return [return_number(c) for c in re.split(r'(\d+)', text) ]
            

#find ORF, return translation
def translate(seq):
	startcodon = re.compile("ATG")
	stopcodon = re.compile("TGA|TAA|TAG")
	ORF = []
	for match_start in startcodon.finditer(seq):
		remaining = seq[match_start.start():]
		if stopcodon.search(remaining):
			for match_stop in stopcodon.finditer(remaining):
				frame = remaining[0:match_stop.end()]
				if len(frame) % 3 == 0:
					ORF.append(frame)
					break
	sorted_ORF = sorted(ORF, key = len)
	final_ORF = Seq(sorted_ORF[len(sorted_ORF)-1])
	protein = str(final_ORF.translate())
	return protein 

#use dendropy to calculate pairwise distances, output a dictionary with all distances 
def distances(input):
    the_dict = {}
    tree = dendropy.Tree.get_from_path(input, "newick", preserve_underscores=True)
    pdm = treemeasure.PatristicDistanceMatrix(tree)
    for i, t1 in enumerate(tree.taxon_namespace):
        tax1 = str(t1)[1:len(str(t1))-1]
        the_dict[str(tax1)] = {}
        for t2 in tree.taxon_namespace[i+1:]:
            tax2 = str(t2)[1:len(str(t2))-1]
            the_dict[str(tax1)][str(tax2)] = round(float(pdm(t1, t2)),5)
    return the_dict

#rebuild a fasta file starting from a defined dictionary, save as defined output
def rebuild_fasta(input_dict, output_fa):
    filename = str(output_fa) + ".fa"
    fileout = open(filename, "w")
    for key in sorted_keys:
        start = 0
        end = 70
        fileout.write (">"+key+"\n")
        seq = input_dict[key]
        while end <= len(seq):
            fileout.write(seq[start:end]+"\n")
            start = end 
            end = end + 70
        fileout.write(seq[start:len(seq)]+"\n")
    fileout.close()
    return filename

###functions to calculate scores for protein alignment

def getScore(a,b,S):
    if (a,b)  in S:
        return S[a,b]
    elif (b,a) in S:
        return S[b,a] 
    else:
        return int(0)
        
def calculate_score(pep1, pep2, matrix): 
    i=0 
    list = []
    while i<len(pep1): 
        Score = getScore(pep1[i],pep2[i], matrix) 
        list.append(Score) 
        i = i+1 
    return list
 
#converts to phylip for phyml
def convert_to_phylip(filename):
    convert_to_phy = AlignIO.parse(filename,"clustal")
    outfile = str(filename) + ".phy"
    AlignIO.write(convert_to_phy, open(outfile,"w"),"phylip")

#extrapolates distance from distance matrices: since matrices of pw distances are not uniform, 
#this search for distance from a to b, and if not found, returns distance from b to a
def extrapolate_distance(dictionary, tx1, tx2):
    if tx1 in dictionary:
        if tx2 in dictionary[tx1]:
            #print ("primary taxon " + str(tx1) + " and secondary taxon " + str(tx2) + " were found in the right order in dictionary ")
            result = dictionary[taxon1][taxon2]
        elif tx2 in dictionary:
            if tx1 in dictionary[taxon2]:
                #print ("primary taxon " + str(tx1) + " was not found in primary keys but in secondary, so results were computed anyway")
                result = dictionary[taxon2][taxon1]
    else:
        result = 0
    return float(result)


#function to plot x and y or several x and one Y 
def my_plot(title, y, x1, x1lab, x2=None, x2lab=None, x3=None, x3lab=None, x4=None, x4lab=None):

    plt.scatter(x1, y, s = 5, marker = '*', label = x1lab)
    fit = np.polyfit(x1, y, 1)
    p = np.poly1d(fit)
    plt.plot(x1, p(x1), "r--", c = 'b')
    
    if x2 != None and x2lab != None:
        plt.scatter(x2, y, s = 5, marker = 'x', label = x2lab)
        fit = np.polyfit(x2, y, 1)
        p = np.poly1d(fit)
        plt.plot(x2, p(x2), "r:", c = 'y')

    if x3 != None and x3lab != None:
        plt.scatter(x3, y, s = 5, marker = '^', label = x3lab)
        fit = np.polyfit(x3, y, 1)
        p = np.poly1d(fit)
        plt.plot(x3, p(x3), "r-.", c = 'g')
    
    if x4 != None and x4lab != None:
        plt.scatter(x4, y, s = 5, marker = "D", label = x4lab)
        fit = np.polyfit(x4, y, 1)
        p = np.poly1d(fit)
        plt.plot(x4, p(x4), "r-", c = 'r')
    
    plt.legend()
    plt.savefig(title)
    plt.close()
######################### THIS IS WHERE THIS STARTS ###########################

#input statement, define input file 
InputDNA = sys.argv[1]
    
InputDNA_noext = InputDNA[0:len(InputDNA)-3] ##change this to regex 
DNAdict = {}


##locate intro statement here, interface and such



#parse DNA input, generate indexed dictionary with trimmed keys -- WORKS
with open(InputDNA) as f:
    for line in f:
        if len(line) > 1:
            if '>' in line:
                key = line.strip()[1:11]
                DNAdict[key] = ""
            else:
                DNAdict[key] += line.strip() 

#a list to sort the fasta files
sorted_keys = []
for key in DNAdict:
	sorted_keys.append(str(key))
sorted_keys.sort(key=natural_keys)
#print sorted_keys

## translate files, store them in a dictionary. Calculate protein structure info as well
AAdict = {}
PROT_info_dict = {}
for key in DNAdict:
	seq_str = str(DNAdict[key])
	PROT_info_dict[key] = info(seq_str)
	AAdict[key] = translate(seq_str)


## rebuild fasta file with proteins
output_protein_name = str(InputDNA_noext) + "_trans"
AA_fasta_ref = rebuild_fasta(AAdict, output_protein_name)

##run ClustalW with DNA

clustalw_NUC = ClustalwCommandline("clustalw2",infile=InputDNA, type='dna', quiet=True)
os.system(str(clustalw_NUC))


##run ClustalW with Protein

clustalw_AA = ClustalwCommandline("clustalw2",infile=AA_fasta_ref, type='PROTEIN', quiet=True)
os.system(str(clustalw_AA))

alnfile_DNA = str(InputDNA_noext) + ".aln"
alnfile_AA = str(output_protein_name) + ".aln"

###make dictionary to store DNA alignments
alignment_dict_DNA = {}
alignment = list(AlignIO.parse(alnfile_DNA,"clustal"))
for i in range(0, len(alignment[0])):
    id = str(alignment[0][i].id)[0:10]
    seq = str(alignment[0][i]._seq)
    alignment_dict_DNA[id] = seq

###make dictionary to store protein alignments
alignment_dict_AA = {}
alignmentAA = list(AlignIO.parse(alnfile_AA,"clustal"))
for i in range(0, len(alignmentAA[0])):
    id = str(alignmentAA[0][i].id)[0:10]
    seq = str(alignmentAA[0][i]._seq)
    alignment_dict_AA[id] = seq

#compute transitions and transversions, gaps, total score. 
#store in dictionary for each sequence pair

tr_tv_dict = {}
for key in alignment_dict_DNA:
    taxon1 = str(key) 
    tr_tv_dict[taxon1] = {}
    for key2 in alignment_dict_DNA:
        taxon2 = str(key2)
        transitions = 0 
        transversions = 0
        gaps = 0
        matches = 0
        for i in range(0, len(alignment_dict_DNA[key])):
            pair_clean = str(alignment_dict_DNA[key][i]) + str(alignment_dict_DNA[taxon2][i])
            if pair_clean == "CG" or pair_clean == "GC" or pair_clean == "AT" or pair_clean == "TA":
                transversions = transversions + 1 
            elif pair_clean == "AC" or pair_clean == "CA" or pair_clean == "TG" or pair_clean == "GT":
                transitions = transitions + 1
            elif pair_clean == "A-" or pair_clean == "-A"or pair_clean == "C-" or pair_clean == "-C" or pair_clean == "G-" or pair_clean == "-G" or pair_clean == "T-" or pair_clean == "-T":
                gaps = gaps + 1 
            elif pair_clean == "AA" or pair_clean == "GG" or pair_clean == "CC" or pair_clean == "TT":
                matches = matches + 1
        tot_score = (2 * matches) + ((-2) * transitions) + ((-3) * gaps) + ((-4) * transversions) 
        total = transitions + transversions
        if total >  0:
            transitions_of_total = float(transitions) / total
            transversions_of_total = float(transversions) / total
        else: 
            transitions_of_total = 0
            transversions_of_total = 0
        
        tr_tv_dict[taxon1][taxon2] = [transitions, round(transitions_of_total,4), transversions, round(transversions_of_total,4), gaps, matches, tot_score]

###compute BLOSUM62 and PAM60 scores based on matrix 

AA_score_dict = {}
for key1 in alignment_dict_AA:
    AA_score_dict[key1] = {}
    for key2 in alignment_dict_AA:
        blo = MatrixInfo.blosum62
        pam = MatrixInfo.pam60
        prot1 = str(alignment_dict_AA[key1]) 
        prot2 = str(alignment_dict_AA[key2]) 
        
        blosum_list = calculate_score(prot1,prot2,blo)
        pam_list = calculate_score(prot1,prot2,pam)
        finalblo = sum(blosum_list)
        finalpam = sum(pam_list)
        AA_score_dict[key1][key2] = [finalblo, finalpam]

#convert aln files to phylip
convert_to_phylip(alnfile_DNA)
convert_to_phylip(alnfile_AA)

#run phyml on DNA
phyfile_DNA = str(alnfile_DNA) + ".phy"

Phyml_runK = PhymlCommandline(input=phyfile_DNA, model='K80')
os.system(str(Phyml_runK)) 
Files = glob.glob('*tree.txt') 
for file in Files: 
    os.rename(file, 'arbolK80')

Phyml_runG = PhymlCommandline(input=phyfile_DNA, model='GTR')
os.system(str(Phyml_runG))
Files = glob.glob('*tree.txt')
for file in Files:
    os.rename(file,'arbolGTR')

Phyml_runH = PhymlCommandline(input=phyfile_DNA, model='HKY85')
os.system(str(Phyml_runH))
Files = glob.glob('*tree.txt')
for file in Files:
    os.rename(file,'arbolHKY85')

Phyml_runJC = PhymlCommandline(input=phyfile_DNA, model='JC69')
os.system(str(Phyml_runJC))
Files = glob.glob('*tree.txt')
for file in Files:
    os.rename(file,'arbolJC69')

#run phyml on AA

phyfile_AA = str(alnfile_AA) + ".phy"

Phyml_runBLO = PhymlCommandline(input=phyfile_AA, datatype="aa", model='Blosum62')
os.system(str(Phyml_runBLO)) 
Files = glob.glob('*tree.txt') 
for file in Files: 
    os.rename(file, 'arbolBLOSUM')

Phyml_runDH = PhymlCommandline(input=phyfile_AA, datatype="aa", model='Dayhoff')
os.system(str(Phyml_runDH)) 
Files = glob.glob('*tree.txt') 
for file in Files: 
    os.rename(file, 'arbolDayhoff')

##compute pairwise distances, store in dictionary

K80_dict = distances("arbolK80")		
HKY85_dict = distances("arbolHKY85")		
JC69_dict = distances("arbolJC69")
GTR_dict = distances("arbolGTR")
BLOSUM_dict = distances("arbolBLOSUM")
DH_dict	= distances("arbolDayhoff")

#use dictionary to loop through all results, print a result table that looks like
#t1 t2 trans trans_of_total trvers tvers_of_total k80 hky85 jc69 gtr blosum dh deltamw deltapi deltacharge

tr_of_total_list = []
tv_of_total_list = []
nucleotide_score_list = []
BLOSUM62_score_list = []
PAM60_score_list = []
d_K80_list = []
d_HKY85_list = []
d_JC69_list = []
d_GTR_list = []
d_BLOSUM_list = []
d_DH_list = []
delta_mw_list = []
delta_pI_list = []
delta_charge_list = []                

handle = open("results.txt", "w")

if len(sys.argv) == 3:
    reference = int(sys.argv[2]) - 1
    taxon1 = str(sorted_keys[reference])
    for taxon2 in sorted_keys:
        if str(taxon1) != str(taxon2):
            #no if statement because this dictionary contains all pairs, even with itself, even redundant
            tr_of_total = tr_tv_dict[taxon1][taxon2][1]       
            tv_of_total = tr_tv_dict[taxon1][taxon2][3]
            nucleotide_score = tr_tv_dict[taxon1][taxon2][6]
            BLOSUM62_score = AA_score_dict[taxon1][taxon2][0]
            PAM60_score = AA_score_dict[taxon1][taxon2][1]
            d_K80 = extrapolate_distance(K80_dict, taxon1, taxon2)
            d_HKY85 = extrapolate_distance(HKY85_dict, taxon1, taxon2)
            d_JC69 = extrapolate_distance(JC69_dict, taxon1, taxon2)
            d_GTR = extrapolate_distance(GTR_dict, taxon1, taxon2)
            d_BLOSUM = extrapolate_distance(BLOSUM_dict, taxon1, taxon2)
            d_DH = extrapolate_distance(DH_dict, taxon1, taxon2)
            
            #no if statement because this puts together stuff using single protein characteristics, hence no pairs in the dictionary         
            delta_mw = round(float(PROT_info_dict[taxon1][0] - PROT_info_dict[taxon2][0]),4)
            delta_pI = round(float(PROT_info_dict[taxon1][1] - PROT_info_dict[taxon2][1]),4)
            delta_charge = round(float(PROT_info_dict[taxon1][2] - PROT_info_dict[taxon2][2]),4)
            
            tr_of_total_list.append(tr_of_total)
            tv_of_total_list.append(tv_of_total)
            nucleotide_score_list.append(nucleotide_score)
            BLOSUM62_score_list.append(BLOSUM62_score)
            PAM60_score_list.append(PAM60_score)
            d_K80_list.append(d_K80)
            d_HKY85_list.append(d_HKY85)
            d_JC69_list.append(d_JC69)
            d_GTR_list.append(d_GTR)
            d_BLOSUM_list.append(d_BLOSUM)
            d_DH_list.append(d_DH)
            delta_mw_list.append(delta_mw)
            delta_pI_list.append(delta_pI)
            delta_charge_list.append(delta_charge)
                            
            print >> handle, taxon1, taxon2, tr_of_total, tv_of_total, nucleotide_score, BLOSUM62_score, PAM60_score, d_K80, d_HKY85, d_JC69, d_GTR, d_BLOSUM, d_DH, delta_mw, delta_pI, delta_charge
    handle.close()
            
elif len(sys.argv) == 2: 
    for t1 in sorted_keys:
        #print ("taxon 1 is " + str(t1))
        for t2 in K80_dict[t1]:
            taxon1 = t1
            taxon2 = t2
            #print ("taxon 2 is " + str(t2))
            tr_of_total = tr_tv_dict[taxon1][taxon2][1]
            tv_of_total = tr_tv_dict[taxon1][taxon2][3]
            nucleotide_score = tr_tv_dict[taxon1][taxon2][6]
            BLOSUM62_score = AA_score_dict[taxon1][taxon2][0]
            PAM60_score = AA_score_dict[taxon1][taxon2][1]
            d_K80 = extrapolate_distance(K80_dict, taxon1, taxon2)
            d_HKY85 = extrapolate_distance(HKY85_dict, taxon1, taxon2)
            d_JC69 = extrapolate_distance(JC69_dict, taxon1, taxon2)
            d_GTR = extrapolate_distance(GTR_dict, taxon1, taxon2)
            d_BLOSUM = extrapolate_distance(BLOSUM_dict, taxon1, taxon2)
            d_DH = extrapolate_distance(DH_dict, taxon1, taxon2)
                    
            #no if statement because this puts together stuff using single protein characteristics, hence no pairs in the dictionary
            delta_mw = round(float(PROT_info_dict[taxon1][0] - PROT_info_dict[taxon2][0]),4)
            delta_pI = round(float(PROT_info_dict[taxon1][1] - PROT_info_dict[taxon2][1]),4)
            delta_charge = round(float(PROT_info_dict[taxon1][2] - PROT_info_dict[taxon2][2]),4)
                                                                                                                                                                                    
            tr_of_total_list.append(tr_of_total)
            tv_of_total_list.append(tv_of_total)
            nucleotide_score_list.append(nucleotide_score)
            BLOSUM62_score_list.append(BLOSUM62_score)
            PAM60_score_list.append(PAM60_score)
            d_K80_list.append(d_K80)
            d_HKY85_list.append(d_HKY85)
            d_JC69_list.append(d_JC69)
            d_GTR_list.append(d_GTR)
            d_BLOSUM_list.append(d_BLOSUM)
            d_DH_list.append(d_DH)
            delta_mw_list.append(delta_mw)
            delta_pI_list.append(delta_pI)
            delta_charge_list.append(delta_charge)

            print >> handle, taxon1, taxon2, tr_of_total, tv_of_total, nucleotide_score, BLOSUM62_score, PAM60_score, d_K80, d_HKY85, d_JC69, d_GTR, d_BLOSUM, d_DH, delta_mw, delta_pI, delta_charge
    handle.close()


my_plot("distances_NUCvsDelta_MW.pdf", delta_mw_list, d_K80_list, "K80", d_HKY85_list,"HKY85", d_JC69_list,"JC69", d_GTR_list, "GTR")
my_plot("distances_NUCvsDelta_pI.pdf", delta_pI_list, d_K80_list, "K80", d_HKY85_list,"HKY85", d_JC69_list,"JC69", d_GTR_list, "GTR")

my_plot("distances_PROTvsDelta_MW.pdf", delta_mw_list, d_BLOSUM_list, "BLOSUM", d_DH_list,"Dayhoff")
my_plot("distances_PROTvsDelta_pI.pdf", delta_pI_list, d_BLOSUM_list, "BLOSUM", d_DH_list,"Dayhoff")

my_plot("NUCscorevsDelta_MW.pdf", delta_mw_list, nucleotide_score_list, "Nucleotide Score")
my_plot("NUCscorevsDelta_pI.pdf", delta_pI_list, nucleotide_score_list, "Nucleotide Score")

my_plot("PROTscorevsDelta_MW.pdf", delta_mw_list, BLOSUM62_score_list, "BLOSUM62", PAM60_score_list,"PAM60")
my_plot("PROTscorevsDelta_pI.pdf", delta_pI_list, BLOSUM62_score_list, "BLOSUM62", PAM60_score_list,"PAM60") 

my_plot("transitions_of_totalvsDelta_MW.pdf", delta_mw_list, tr_of_total_list, "transitions")
my_plot("transitions_of_totalvsDelta_pI.pdf", delta_pI_list, tr_of_total_list, "transitions")

my_plot("transversions_of_totalvsDelta_MW.pdf", delta_mw_list, tv_of_total_list, "transversions")
my_plot("transversions_of_totalvsDelta_pI.pdf", delta_pI_list, tv_of_total_list, "transversions")
