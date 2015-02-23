import re
import sys
import pickle

##########################################################################################################################
# INPUT FILES
##########################################################################################################################
print "##gff-version 3"
#############################################
# Step I: Finding Scaffold lengths 
#############################################
# Length of Sccfolds in Salmon 3.6 Chr: Updated 25th Aug 2014 
# This file contain length of scaffolds in Salmon 3.6 Chr. Its an output from "samtools faidx XXXX.fasta"  
# Updated 25th of Aug 2014
SL1 = open("Salmon_3p6_Chr_Data/CIGENE_3.6chr-with_unmapped_genes.fasta.fai","r")
SL = {}
for i in SL1:
    temp = i.split()
    SL[temp[0]] = int(temp[1])

###########################################
# Step 2: New names for gene and range.
###########################################

# Only the longest Locus with support
# open("CIGSSA_Longest_Support.annot.tab","r")

#  ALL alternative Transcripts in each locus
# open("CIGSSA_All.annot.tab","r")


# New blocks of code incorporated for new file format from Simen on 25th Aug 2014. 
#FD1_new = open("CIGSSA_Longest_All_Complete.annot.tab","r")

# Updated File : 120914
FD1_new = open("Salmon_3p6_Chr_Data/CIGSSA_All.annot-withsupport.tab","r")
FD_new = FD1_new.readlines()

GeneNames_New = {}
GeneRange_New = {}
GeneModels = {}
GeneModels_Gene = {}
Gene_Scaffold = {}
Trans_Range = {}
# Gene/transcript orientation for Reference strand
GOrin = {}

# Gene/transcript orientation for Trans Decoder strand
GOrin_TD = {}


for i in xrange(1,len(FD_new)):
    temp = FD_new[i].split()
    GeneNames_New[temp[0]] = temp[14]
    #GeneRange_New[temp[0]] = str(temp[-2]) + "-" + str(temp[-1].strip())
    GeneRange_New[temp[0]] = str(temp[-4]) + "-" + str(temp[-3].strip())
    Trans_Range[temp[0]] = str(temp[-2]) + "-" + str(temp[-1].strip())
    GeneModels[temp[5]] = temp[5]
    ID = temp[5].replace("m","g")
    GeneModels_Gene[temp[5]] = ID
    Gene_Scaffold[temp[0]] = temp[1]
    
    # Storing the Gene/Transcript "REFERENCE" orientation 
    GOrin[temp[0]] = temp[2]

    # Storing the Gene/Transcript "Trans Decoder" orientation 
    GOrin_TD[temp[0]] = temp[4]

##########################################
# Step 3: Exon information 
##########################################
# File was updated on 25th Aug 2014
MG1 = open("Salmon_3p6_Chr_Data/merged.gtf","r")
#MG1 = open("Test_Merged.gtf","r")
MG = MG1.readlines()

# Gene, mRNA, UTR, CDS information:
# Input File Updated 25 Aug 2014
GM1 = open("Salmon_3p6_Chr_Data/TD_concatenated.gff","r")
#GM1 = open("Test_TD.gff","r")
GM = GM1.readlines()


####################################################
# Step 5: List of scaffold names in Salmon 3.6 Chr 
###################################################
# Read the scaffodls available in Salmon 3.6 Chr assembly
# Updated 25 Aug 2014
All_Scaf_List1 = open("Salmon_3p6_Chr_Data/Salmon_3p6_Chr_Scaf_Names_250814.txt","r")

# Store the results in a dict
All_Scaf_List = {}
for i in All_Scaf_List1:
    ID = i.replace(">","")
    All_Scaf_List[ID.strip()] = ID.strip()

for i in All_Scaf_List:
    ID = i.strip()
    txt_1 = ID + "\t" + "Salmon" + "\t" + "sequence_assembly" + "\t" + str(1) + "\t" + str(SL[ID]) + "\t" + "." + "\t" + "+" + "\t" + "." + "\t" + "ID=" + ID + ";" + "Name=" + i + "\n"
    print txt_1.strip()
    
#############################################################################################################################
# Precited gene names from FABIAN by end of May 2014
#############################################################################################################################
#GN1 = open("/mnt/users/jeevka/Sally_Gene_Names_Function/Sally_Gene_Names_061014.txt","r")
GN1 = open("/mnt/users/fabig/Ssa_transcriptome/cigene3.6_chrom/protein_db/Final_blast_uniprot.csv","r")
GN = GN1.readlines()

GType = {}
GSymbol = {}
for i in xrange(1,len(GN)):
    temp = GN[i].split(",")
    ID = temp[1].replace("\"","")
    symbol = temp[3].replace("\"","")
    func = temp[2].replace("\"","")
    
    # Storing the gene symbol
    GSymbol[ID] = symbol
    
    # Storing predicted Gene function
    """
    txt = ""
    for j in range(1,len(temp)-1):
        txt += temp[j] + " "
    """
    
    GType[ID] = func

# Predicted GO term for the transcripts
GOT1 = open("/mnt/users/fabig/Ssa_transcriptome/cigene3.6_chrom/data/cig36chrom_b2go_annot.txt","r")
GOT = GOT1.readlines()
GOTerm = {}
for i in xrange(1,len(GOT)):
    ID = GOT[i].split()[0]
    temp1 = GOT[i].split("GO:")
    temp2 = temp1[1].split()
    txt = ""
    for j in xrange(1,len(temp2)):
        txt += temp2[j] + " "
    
    GOTerm[ID] = txt.strip()

##########################################################################################################################
# Predicted Gene names: Torfinn on 06-06-14
##########################################################################################################################
"""
# HAVE TO WAIT FOR THE UPDATED NEW FILE FROM TORFINN/SIMEN
GN1 = open("/mnt/backup2/users/aquagenome/Ssa/Synteny/CIGENE-3.6-unmasked/ANOT/From-Zebrafish/ssaCIG_homolyevidence_all.transcripts_GBROWSEinfo-onlyT1s-with-annotation.tab","r")
GN = GN1.readlines()

GType = {}; GSymbol = {}
for i in GN:
    temp = i.split()
    GSymbol[temp[1]] = temp[-1]
    function = temp[-2].replace(";","-")
    GType[temp[1]] = function

for i in GType:
    print i,"\t",GType[i],"\t",GSymbol[i]
sys.exit()
"""

###########################################################################################################################
#  SUB - PROGRAMS 
###########################################################################################################################
def store_data_1(TRANS,TID,KEY,START,END):
    data = str(START) + "-" + str(END)
    TRANS[TID][KEY] = data

    return TRANS

def store_data_2(TRANS,TID,KEY,START,END):
    data = str(START) + "-" + str(END)
    TRANS[TID] = {KEY: data}

    return TRANS

def separate_range(RANGE):
    temp = RANGE.split("-")
    
    return int(temp[0]),int(temp[1])

###########################################################################################################################
# Decide whether Cufflinks predicted Exons fall within CDS (Transdecode) range
# UTR ranges (start and end) and exons ranges (R1 and R2) and whether UTRs are present or not?
# if UTR5 == 1, 5 UTR is present. If UTR% == 0, 5 UTR is not present
###########################################################################################################################
def decide_on_CDS(start,end,R1,R2,UTR5,UTR3,strand,GeneName):
    # If its "+" strand
    # Start: end of 5UTR.
    # end: begnning of 3UTR
    
    # If its "-" strand
    # Start: end of 3UTR.
    # end: begnning of 5UTR
    
    # Indicator wherther CDS is ok or not
    CDS = 0
    
    # Case 1: Exon boundaries falls within UTR ranges
    if int(R1) > start and int(R2) < end:
        #print "Case I:",GeneName,"\t",R1,"\t",R2,"\t",start,"\t",end        
        cds_start = R1; cds_end = R2
        CDS = 1
    
    # Case II: Exon boundaries are equal to UTR boundaries    
    if int(R1) == start and int(R2) == end:
        #print "Case II:",GeneName,"\t",R1,"\t",R2,"\t",start,"\t",end
        if UTR5 == 1:
            cds_start = R1 + 1
        else:
            cds_start = R1
        
        if UTR3 == 1:
            cds_end = R2 - 1
        else:
            cds_end = R2
            
        CDS = 1
    
    # Case III: Exon start is less than UTR start and Exon end is bigger than UTR start 
    if int(R1) < start and int(R2) > start and int(R2) < end:
        #print "Case III:",GeneName,"\t",R1,"\t",R2,"\t",start,"\t",end
        if strand == "+":
            if UTR5 == 1 and UTR3 == 1:
                cds_start = start + 1
                cds_end = R2
            elif UTR5 == 1 and UTR3 == 0:
                cds_start = start + 1
                cds_end = R2
            elif UTR5 == 0 and UTR3 == 1:
                cds_start = start
                cds_end = R2
            else:
                cds_start = start
                cds_end = R2            
        else:
            if UTR5 == 1 and UTR3 == 1:
                cds_start = start + 1
                cds_end = R2
            elif UTR5 == 1 and UTR3 == 0:
                cds_start = start
                cds_end = R2
            elif UTR5 == 0 and UTR3 == 1:
                cds_start = start + 1
                cds_end = R2
            else:
                cds_start = start
                cds_end = R2
                
        CDS = 1        
    
    # Case IV: Exon start is bigger than UTR start and exon end is also bigger than UTR end
    if int(R1) > start and int(R1) < end and int(R2) > end:
        #print "Case IV:",GeneName,"\t",R1,"\t",R2,"\t",start,"\t",end
        if strand == "+":
            if UTR5 == 1 and UTR3 == 1:
                cds_start = R1
                cds_end = end - 1
            elif UTR5 == 1 and UTR3 == 0:
                cds_start = R1
                cds_end = end
            elif UTR5 == 0 and UTR3 == 1:
                cds_start = R1
                cds_end = end - 1
            else:
                cds_start = R1
                cds_end = end                
        else:
            if UTR5 == 1 and UTR3 == 1:
                cds_start = R1
                cds_end = end - 1
            elif UTR5 == 1 and UTR3 == 0:
                cds_start = R1
                cds_end = end - 1
            elif UTR5 == 0 and UTR3 == 1:
                cds_start = R1
                cds_end = end
            else:
                cds_start = R1
                cds_end = end        
        CDS = 1
    
    # Case V: Exon start is less than UTR start and exon end is equal to UTR end
    if int(R1) < start  and int(R2) == end:
        #print "Case V",GeneName,"\t",R1,"\t",R2,"\t",start,"\t",end
        if strand == "+":
            if UTR5 == 1 and UTR3 == 1:
                cds_start = start + 1
                cds_end = end - 1
            elif UTR5 == 1 and UTR3 == 0:
                cds_start = start + 1
                cds_end = end
            elif UTR5 == 0 and UTR3 == 1:
                cds_start = start
                cds_end = end - 1
            else:
                cds_start = start
                cds_end = end                
        else:
            if UTR5 == 1 and UTR3 == 1:
                cds_start = start + 1
                cds_end = end - 1
            elif UTR5 == 1 and UTR3 == 0:
                cds_start = start
                cds_end = end - 1
            elif UTR5 == 0 and UTR3 == 1:
                cds_start = start + 1
                cds_end = end
            else:
                cds_start = start
                cds_end = end        
        CDS = 1

    # Case VI: Exon start is equal to UTR start and exon end is bigger than UTR end
    if int(R1) == start  and int(R2) > end:
        #print "Case VI",GeneName,"\t",R1,"\t",R2,"\t",start,"\t",end
        if strand == "+":
            if UTR5 == 1 and UTR3 == 1:
                cds_start = start + 1
                cds_end = end - 1
            elif UTR5 == 1 and UTR3 == 0:
                cds_start = start + 1
                cds_end = end
            elif UTR5 == 0 and UTR3 == 1:
                cds_start = start
                cds_end = end - 1
            else:
                cds_start = start
                cds_end = end                
        else:
            if UTR5 == 1 and UTR3 == 1:
                cds_start = start + 1
                cds_end = end - 1
            elif UTR5 == 1 and UTR3 == 0:
                cds_start = start
                cds_end = end - 1
            elif UTR5 == 0 and UTR3 == 1:
                cds_start = start + 1
                cds_end = end
            else:
                cds_start = start
                cds_end = end
        
        CDS = 1

    # Case VII: Exon start is less than UTR start and exon end is bigger than UTR end
    if int(R1) < start  and int(R2) > end:
        #print "Case VII:",GeneName,"\t",R1,"\t",R2,"\t",start,"\t",end
        if strand == "+":
            if UTR5 == 1 and UTR3 == 1:
                cds_start = start + 1
                cds_end = end - 1
            elif UTR5 == 1 and UTR3 == 0:
                cds_start = start + 1 
                cds_end = end
            elif UTR5 == 0 and UTR3 == 1:
                cds_start = start
                cds_end = end - 1
            else:
                cds_start = start
                cds_end = end                
        else:
            if UTR5 == 1 and UTR3 == 1:
                cds_start = start + 1 
                cds_end = end - 1
            elif UTR5 == 1 and UTR3 == 0:
                cds_start = start
                cds_end = end - 1
            elif UTR5 == 0 and UTR3 == 1:
                cds_start = start + 1 
                cds_end = end
            else:
                cds_start = start
                cds_end = end        

        CDS = 1

    # Case IV: Exon start is bigger than UTR start and exon end is also bigger than UTR end
    if int(R1) > start and int(R1) < end and int(R2) == end:
        #print "Case VIII:",GeneName,"\t",R1,"\t",R2,"\t",start,"\t",end
        if strand == "+":
            if UTR5 == 1 and UTR3 == 1:
                cds_start = R1
                cds_end = end - 1
            elif UTR5 == 1 and UTR3 == 0:
                cds_start = R1
                cds_end = end
            elif UTR5 == 0 and UTR3 == 1:
                cds_start = R1
                cds_end = end - 1 
            else:
                cds_start = R1
                cds_end = end                
        else:
            if UTR5 == 1 and UTR3 == 1:
                cds_start = R1
                cds_end = end - 1
            elif UTR5 == 1 and UTR3 == 0:
                cds_start = R1
                cds_end = end - 1 
            elif UTR5 == 0 and UTR3 == 1:
                cds_start = R1
                cds_end = end
            else:
                cds_start = R1
                cds_end = end        
        CDS = 1        

    # Case IX: Exon start is equal to UTR start and exon end is bigger than UTR end
    if int(R1) == start  and int(R2) < end:
        #print "Case IX",GeneName,"\t",R1,"\t",R2,"\t",start,"\t",end
        if strand == "+":
            if UTR5 == 1 and UTR3 == 1:
                cds_start = start + 1
                cds_end = R2
            elif UTR5 == 1 and UTR3 == 0:
                cds_start = start + 1
                cds_end = R2
            elif UTR5 == 0 and UTR3 == 1:
                cds_start = start
                cds_end = R2
            else:
                cds_start = start
                cds_end = R2                
        else:
            if UTR5 == 1 and UTR3 == 1:
                cds_start = start + 1
                cds_end = R2
            elif UTR5 == 1 and UTR3 == 0:
                cds_start = start
                cds_end = R2
            elif UTR5 == 0 and UTR3 == 1:
                cds_start = start + 1 
                cds_end = R2
            else:
                cds_start = start
                cds_end = R2
        
        CDS = 1

    return cds_start,cds_end,CDS

def store_gene_ranges(Longest_Gene_Range_Start,Longest_Gene_Range_End,trans,gstart,gend):
    
    # If the gene is in positive strand
    if gstart < gend:                 
        if Longest_Gene_Range_Start.has_key(trans):
            if Longest_Gene_Range_Start[trans] > gstart:
                Longest_Gene_Range_Start[trans] = gstart
                
            if Longest_Gene_Range_End[trans] < gend:
                Longest_Gene_Range_End[trans] = gend
        else:
            Longest_Gene_Range_Start[trans] = gstart
            Longest_Gene_Range_End[trans] = gend
            
    # If the gene is in Negative strand
    else:
        if Longest_Gene_Range_Start.has_key(trans):
            if Longest_Gene_Range_Start[trans] < gstart:
                Longest_Gene_Range_Start[trans] = gstart
                
            if Longest_Gene_Range_End[trans] > gend:
                Longest_Gene_Range_End[trans] = gend
        else:
            Longest_Gene_Range_Start[trans] = gstart
            Longest_Gene_Range_End[trans] = gend        
    
    return Longest_Gene_Range_Start,Longest_Gene_Range_End
    
#######################################################################################################
# This function is needed because the Trans decoder considers
# the UTR regions as coding regions.
# E.g. If 5 UTR is 5000bps from TD results, its only considers 5000 coding base pairs
# Which can span for more than one or two exons.
# We should skip the introns in those regions.
#######################################################################################################
#######################################################################################################
# UTR Regions can span for more than one exon. This funtion will take care of that situations for 5 UTR
#######################################################################################################
def calculate_5_UTR_regions(start,end,EXONS,gstart,gend,GOrin,GOrin_TD,trans,GeneNames_New,temp_exon):
    
    # Decide the strand 
    # Gene orientation
    if GOrin[trans] == ".":
        Gene_ori = GOrin_TD[trans]
        GO = GOrin_TD[trans]
    else:
        Gene_ori = GOrin[trans]
        GO = GOrin[trans]
        
    txt_5UTR = "" 
    
    # Actual UTR Size
    length_of_UTR = end - start + 1
    
    # If its a "+" strand or ". (no info)"
    if Gene_ori == "+" or Gene_ori == ".":
        # New UTR Size 
        UTR_Size = 0
        for i in EXONS:
            temp = EXONS[i].split("-")
            ex_size = int(temp[1]) - int(temp[0]) + 1
            
            if GOrin[trans] == ".":
                GO = GOrin_TD[trans]
            else:
                GO = GOrin[trans]            
            
            if (UTR_Size + ex_size) < length_of_UTR:
                txt_5UTR += scaf + "\t" + "Salmon" + "\t" + "five_prime_UTR" + "\t" + str(temp[0]) + "\t" + str(temp[1]) + "\t" + str(temp_exon[5]) + "\t" + GO + "\t" + str(temp_exon[7]) + "\t" + "Parent=" + GeneNames_New[trans] +"\n"
                UTR_Size += ex_size
                
            else:
                temp_UTR_size = UTR_Size + ex_size
                diff = temp_UTR_size - length_of_UTR
                UTR_END = int(temp[1]) - diff
                
                txt_5UTR += scaf + "\t" + "Salmon" + "\t" + "five_prime_UTR" + "\t" + str(temp[0]) + "\t" + str(UTR_END) + "\t" + str(temp_exon[5]) + "\t" + GO + "\t" + str(temp_exon[7]) + "\t" + "Parent=" + GeneNames_New[trans] +"\n"
                break
        
        UTR_START = gstart
    
    # If its "-" strand
    if Gene_ori == "-":    
        # New UTR Size 
        UTR_Size = 0
        for i in xrange(len(EXONS),0,-1):
            temp = EXONS[i].split("-")
            ex_size = int(temp[1]) - int(temp[0]) + 1
            
            if (UTR_Size + ex_size) < length_of_UTR:
                UTR_Size += ex_size
                txt_5UTR += scaf + "\t" + "Salmon" + "\t" + "five_prime_UTR" + "\t" + str(temp[0]) + "\t" + str(temp[1]) + "\t" + str(temp_exon[5]) + "\t" + GO + "\t" + str(temp_exon[7]) + "\t" + "Parent=" + GeneNames_New[trans] +"\n"
                
            else:
                temp_UTR_size = UTR_Size + ex_size
                diff = temp_UTR_size - length_of_UTR
                UTR_START = int(temp[0]) + diff #+ 1
                txt_5UTR += scaf + "\t" + "Salmon" + "\t" + "five_prime_UTR" + "\t" + str(UTR_START) + "\t" + str(temp[1]) + "\t" + str(temp_exon[5]) + "\t" + GO + "\t" + str(temp_exon[7]) + "\t" + "Parent=" + GeneNames_New[trans] +"\n"
                break
            
        UTR_END = gend

    return UTR_START,UTR_END,txt_5UTR

#######################################################################################################
# UTR Regions can span for more than one exon. This funtion will take care of that situations for 3 UTR
#######################################################################################################
def calculate_3_UTR_regions(start,end,EXONS,gstart,gend,GOrin,GOrin_TD,trans,GeneNames_New,temp_exon):
    
    # Decide the strand 
    # Gene orientation
    if GOrin[trans] == ".":
        Gene_ori = GOrin_TD[trans]
        GO = GOrin_TD[trans]
    else:
        Gene_ori = GOrin[trans]
        GO = GOrin[trans]
        
    # Actual UTR Size
    length_of_UTR = end - start + 1
    
    txt_3UTR = ""
    
    # If its a "+" strand or ". (no info)"
    if Gene_ori == "+" or Gene_ori == ".":
        # New UTR Size 
        UTR_Size = 0
        for i in xrange(len(EXONS),0,-1):
            temp = EXONS[i].split("-")
            ex_size = int(temp[1]) - int(temp[0]) + 1
            
            if (UTR_Size + ex_size) < length_of_UTR:
                UTR_Size += ex_size
                txt_3UTR += Gene_Scaffold[trans] + "\t" + "Salmon" + "\t" + "three_prime_UTR" + "\t" + str(temp[0]) + "\t" + str(temp[1]) + "\t" + str(temp_exon[5]) + "\t" + GO + "\t" + str(temp_exon[7]) + "\t" + "Parent=" + GeneNames_New[trans] +"\n"
            else:
                temp_UTR_size = UTR_Size + ex_size
                diff = temp_UTR_size - length_of_UTR
                UTR_START = int(temp[0]) + diff #+ 1
                if UTR_START > int(temp[1]):
                    temp[1] = UTR_START
                    
                txt_3UTR += Gene_Scaffold[trans] + "\t" + "Salmon" + "\t" + "three_prime_UTR" + "\t" + str(UTR_START) + "\t" + str(temp[1]) + "\t" + str(temp_exon[5]) + "\t" + GO + "\t" + str(temp_exon[7]) + "\t" + "Parent=" + GeneNames_New[trans] +"\n"
                break
        
        UTR_END = gend
    
    # If its a "+" strand or ". (no info)"
    if Gene_ori == "-" :
        # New UTR Size 
        UTR_Size = 0
        for i in EXONS:
            temp = EXONS[i].split("-")
            ex_size = int(temp[1]) - int(temp[0]) + 1
            
            if (UTR_Size + ex_size) < length_of_UTR:
                UTR_Size += ex_size
                txt_3UTR += Gene_Scaffold[trans] + "\t" + "Salmon" + "\t" + "three_prime_UTR" + "\t" + str(temp[0]) + "\t" + str(temp[1]) + "\t" + str(temp_exon[5]) + "\t" + GO + "\t" + str(temp_exon[7]) + "\t" + "Parent=" + GeneNames_New[trans] +"\n"
                #print GeneNames_New[trans],"\t",UTR_Size
            else:
                temp_UTR_size = UTR_Size + ex_size
                #print GeneNames_New[trans],"\t",temp_UTR_size
                diff = temp_UTR_size - length_of_UTR
                #print GeneNames_New[trans],"\t",diff
                UTR_END = int(temp[1]) - diff #+ 1 
                txt_3UTR += Gene_Scaffold[trans] + "\t" + "Salmon" + "\t" + "three_prime_UTR" + "\t" + str(temp[0]) + "\t" + str(UTR_END) + "\t" + str(temp_exon[5]) + "\t" + GO + "\t" + str(temp_exon[7]) + "\t" + "Parent=" + GeneNames_New[trans] +"\n"
                break
        
        UTR_START = gstart

    return UTR_START,UTR_END,txt_3UTR


###########################################################################
# In case of 5 UTR missing
# Trans decoder structures are differentfor + and - strands
# when Ref is + and TD is -, things are very different
# This function will take care of that situation for 5 UTR
########################################################################### 
def calculate_5_UTR_regions_in_case_of_missing(TRANS,TD_GO,GO,gstart,gend):
                    
    temp_CDS = TRANS["CDS"].split("-")
    temp_EXON = TRANS["exon"].split("-")

    diff_1 = int(temp_CDS[0]) - int(temp_EXON[0])
    diff_2 = int(temp_EXON[1]) - int(temp_CDS[1])
    
    # Check whether the reference strand and TD strand are same or different
    if GO == "+" and TD_GO == "+":
        UTR5_start = gstart + diff_1
        UTR5_end = gstart  + diff_1
    
    if GO == "-" and TD_GO == "-":
        UTR5_start = gend - diff_2
        UTR5_end = gend - diff_2

    if GO == "+" and TD_GO == "-":
        UTR5_start = gstart + diff_2
        UTR5_end = gstart  + diff_2

    if GO == "-" and TD_GO == "+":
        UTR5_start = gend - diff_1
        UTR5_end = gend - diff_1
                      
    return UTR5_start,UTR5_end

###########################################################################
# In case of 3 UTR missing
# Trans decoder structures are differentfor + and - strands
# when Ref is + and TD is -, things are very different
# This function will take care of that situation for 3 UTR
###########################################################################
def calculate_3_UTR_regions_in_case_of_missing(TRANS,TD_GO,GO,gstart,gend):

    temp_EXON = TRANS["exon"].split("-")
    temp_CDS = TRANS["CDS"].split("-")
                    
    diff_1 = int(temp_CDS[0]) - int(temp_EXON[0])
    diff_2 = int(temp_EXON[1]) - int(temp_CDS[1])

    # Check whether the reference strand and TD strand are same or different
    if GO == "+" and TD_GO == "+":
        UTR3_start = gend - diff_2
        UTR3_end = gend - diff_2
    
    if GO == "-" and TD_GO == "-":
        UTR3_start = gstart + diff_1
        UTR3_end = gstart + diff_1

    if GO == "+" and TD_GO == "-":
        UTR3_start = gend - diff_1
        UTR3_end = gend  - diff_1

    if GO == "-" and TD_GO == "+":
        UTR3_start = gstart + diff_2
        UTR3_end = gstart + diff_2                  
        
    return UTR3_start,UTR3_end

###########################################################################################################################
# OUTPUT FILE
###########################################################################################################################
SGFF = open("Sally_Transcriptome_250814.gff3","w")

# First line of any GFF3 File
SGFF.write("##gff-version 3\n")
#print "##gff-version 3"

###########################################################################################################################
# PART II: Analyze mRNA, CDS, UTR REGION DETAILS
###########################################################################################################################
TRANS = {}

for i in GM:
    temp = i.split()
    
    # ID in FT:  Check whether the ID is in Filtered list
    if len(temp) > 0:
        
        if re.search("CDS",i):
            temp1 = temp[8].split(";")
            ID =  temp1[0].split("=")[1]
            ID1 = ID.split(".")
            ID = ID1[1] + "." + ID1[2]
        else:
            temp1 = temp[8].split(";")
            ID1 =  temp1[0].split("=")[1]
            ID2 =  ID1.split(".")
            ID = ID2[0] + "." + ID2[1]
        
        if re.search("gene",i):
            if GeneModels_Gene.has_key(ID):
                if TRANS.has_key(temp[0]):
                    TRANS = store_data_1(TRANS,temp[0],temp[2],temp[3],temp[4])
                else:
                    TRANS = store_data_2(TRANS,temp[0],temp[2],temp[3],temp[4])
        else:
            if GeneModels.has_key(ID):
                if TRANS.has_key(temp[0]):
                    TRANS = store_data_1(TRANS,temp[0],temp[2],temp[3],temp[4])
                else:
                    TRANS = store_data_2(TRANS,temp[0],temp[2],temp[3],temp[4])

###########################################################################################################################
# PART III: Preprocess the GFF file to get the Gene-Boundary
###########################################################################################################################
temp = MG[0].split()
temp2  = temp[11].replace("\"","")
TID  = temp2.replace(";","")
start = temp[3]

# Gene orientation 
# GOrin = {}

# All the Exons in the gene 
EXONS = {}

# Gene Range
GRange = {}

for i in MG:
    temp = i.split()
    temp2  = temp[11].replace("\"","")
    TID_temp  = temp2.replace(";","")
    
    # Store Exons
    #EXONS[TID_temp] = str(temp[3]) + "-" + str(temp[4])
    
    if EXONS.has_key(TID_temp):
        L = len(EXONS[TID_temp]) +1
        txt = str(temp[3]) + "-" + str(temp[4])
        EXONS[TID_temp][L] = txt
    else:
        txt = str(temp[3]) + "-" + str(temp[4])
        EXONS[TID_temp] = {1: txt}
        
    # Store Gene orientation
    # GOrin[TID_temp] = temp[6]
    
    # Check whether this gene model passed the filter
    if TID != TID_temp:
        data = str(start) + "-" + str(end)
        GRange[TID] = data
        start = temp[3]
        TID = TID_temp
    
    end = temp[4]


GRange[TID] = data

###########################################################################################################################
# PART V: Analyse GTF Files which contains the Exons details
###########################################################################################################################
temp = MG[0].split()
temp2  = temp[11].replace("\"","")
TID  = temp2.replace(";","")

N_Ex = 0
txt_CDS = ""
MM = 0
N_3UTR = 0
N_5UTR = 0

Scaf_List = {}
NM_Scaf = 0

Gene_Lines = {}
Longest_Gene_Range_Start = {}
Longest_Gene_Range_End = {}
Gene_Scaf = {}
Gene_Strand = {}

for i in MG:
    temp = i.split()
    temp2  = temp[11].replace("\"","")
    TID_temp  = temp2.replace(";","")
    scaf = temp[0]
    
    if TRANS.has_key(TID_temp): # and TID_temp in FT:
        # Check whether this gene model passed the filter
        if TID != TID_temp:
                        
            # Printing the CDS/Exon Part
            temp_3 = txt_CDS.split("\t")
            if len(temp_3) > 1:
                print txt_CDS.strip()
            
            SGFF.write(txt_CDS)
            
            ###############################
            # Printing THREE prime UTR
            ###############################            
            try:
                if TRANS[TID].has_key("three_prime_UTR"):
                    UTR3 = 1
                else:
                    UTR3 = 0
                    
                start,end = separate_range(TRANS[TID]["three_prime_UTR"])
                UTR3_start,UTR3_end,txt_3UTR = calculate_3_UTR_regions(start,end,EXONS[trans],gstart,gend,GOrin,GOrin_TD,trans,GeneNames_New,temp)
                #########################################################
                # Irrespective of the strand orientation (+ or -)       #
                # start and end of features are mentioned/considered    #
                # as + oriented to avoid confusion                      #
                #########################################################                
                # Gene orientation
                # Reference and TD strand information is differnet in many cases.                
                #txt_3UTR = Gene_Scaffold[trans] + "\t" + "Salmon" + "\t" + "three_prime_UTR" + "\t" + str(UTR3_start) + "\t" + str(UTR3_end) + "\t" + str(temp[5]) + "\t" + GOrin[trans] + "\t" + str(temp[7]) + "\t" + "Parent=" + GeneNames_New[trans] +"\n"
                print txt_3UTR.strip()        
                
                SGFF.write(txt_3UTR)
            except:
                N_3UTR += 1
                pass
            
            txt_CDS = ""
            
            N_Ex = 0
            TID = TID_temp
          
        if TID == TID_temp:
            N_Ex += 1
            
            # WRITING THE SCAFFOLD INFORMATION AS THE FIRST DATA LINE before the first exon
            if N_Ex == 1:
                    trans = TID
                    
                    ################################
                    # Printing sequence assembly part
                    ################################
                    """
                    try:
                        txt_1 = scaf + "\t" + "Salmon" + "\t" + "sequence_assembly" + "\t" + str(1) + "\t" + str(SL[scaf]) + "\t" + "." + "\t" + "+" + "\t" + "." + "\t" + "ID=" + scaf + ";" + "Name=" + scaf + "\n"
                        print txt_1.strip()
                        
                        Scaf_List[scaf] = scaf
                        
                        SGFF.write(txt_1)        
                    except:
                        NM_Scaf += 1
                        pass
                    """
                    ################################
                    # Printing Gene part
                    ################################
                    gstart,gend = separate_range(GeneRange_New[TID])
                    
                    # store the gene ranges to find the longest ones
                    Longest_Gene_Range_Start,Longest_Gene_Range_End = store_gene_ranges(Longest_Gene_Range_Start,Longest_Gene_Range_End,trans,gstart,gend)
                    Gene_Scaf[trans] = scaf
                    Gene_Strand[trans] = temp[6]
                    
                    if GOrin[trans] == ".":
                        GO = GOrin_TD[trans]
                    else:
                        GO = GOrin[trans]
                    
                    try:
                        GNN = GeneNames_New[trans].split(".")[0]
                        GNN1 = GeneNames_New[trans]
                        txt_gene = scaf + "\t" + "Salmon" + "\t" + "gene" + "\t" + str(gstart) + "\t" + str(gend) + "\t" + str(temp[5]) + "\t" + GO + "\t" + str(temp[7]) + "\t" + "ID=" + GNN + ";Name=" + GNN + ";Alias=" + GSymbol[GNN] + ";Note="+ GType[GNN] + "(predicted)" + "\n"
                        txt_gene_1 = scaf + "\t" + "Salmon" + "\t" + "gene" + "\t" + str(gstart) + "\t" + str(gend) + "\t" + str(temp[5]) + "\t" + GO + "\t" + str(temp[7]) + "\t" + "ID=" + GNN + ";Name=" + GNN + ";Alias=" + GType[GNN] + ";Note="+ GSymbol[GNN] + "(predicted)" + "\n"
                    except:
                        GNN = GeneNames_New[trans].split(".")[0]
                        GNN1 = GeneNames_New[trans]
                        txt_gene = scaf + "\t" + "Salmon" + "\t" + "gene" + "\t" + str(gstart) + "\t" + str(gend) + "\t" + str(temp[5]) + "\t" + GO + "\t" + str(temp[7]) + "\t" + "ID=" + GNN + ";Name=" + GNN + ";Alias=" + "Unknown" + ";Note="+ "Unknown" + "(predicted)" + "\n"
                        txt_gene_1 = scaf + "\t" + "Salmon" + "\t" + "gene" + "\t" + str(gstart) + "\t" + str(gend) + "\t" + str(temp[5]) + "\t" + GO + "\t" + str(temp[7]) + "\t" + "ID=" + GNN + ";Name=" + GNN + "\n"

                    if not Gene_Lines.has_key(GNN):
                        Gene_Lines[GNN] = txt_gene

                    SGFF.write(txt_gene)

                    ################################
                    # Printing mRNA part
                    ################################
                    mRNA_ID = trans + ".1"
                    GNN = GeneNames_New[trans].split(".")[0]
                    GNN1 = GeneNames_New[trans]
                    
                    if GOrin[trans] == ".":
                        GO = GOrin_TD[trans]
                    else:
                        GO = GOrin[trans]
		    try:
	                    txt_mRNA = scaf + "\t" + "Salmon" + "\t" + "mRNA" + "\t" + str(gstart) + "\t" + str(gend) + "\t" + str(temp[5]) + "\t" + GO + "\t" + str(temp[7]) + "\t" + "ID=" + GeneNames_New[trans] + ";Parent=" + GNN + ";Name=" + GeneNames_New[trans] + ";Note=" + GType[GNN] + ";Alias=" + GSymbol[GNN] + ";Ontology_term="+ GOTerm[GNN1] + "\n"
		    except: 
			    txt_mRNA = scaf + "\t" + "Salmon" + "\t" + "mRNA" + "\t" + str(gstart) + "\t" + str(gend) + "\t" + str(temp[5]) + "\t" + GO + "\t" + str(temp[7]) + "\t" + "ID=" + GeneNames_New[trans] + ";Parent=" + GNN + ";Name=" + GeneNames_New[trans] + ";Note=" + "Unknown"  + ";Alias=" + "Unknown" + ";Ontology_term=" + "Unknown" + "\n"
                    print txt_mRNA.strip()

                    SGFF.write(txt_mRNA)

                    ################################
                    # Printing FIVE prime UTR
                    ################################
                    try:
                        start,end = separate_range(TRANS[TID]["five_prime_UTR"])
                        UTR5_start,UTR5_end,txt_5UTR = calculate_5_UTR_regions(start,end,EXONS[trans],gstart,gend,GOrin,GOrin_TD,trans,GeneNames_New,temp)

                        #########################################################
                        # Irrespective of the strand orientation (+ or -)       #
                        # start and end of features are mentioned/considered    #
                        # as + oriented to avoid confusion                      #
                        #########################################################
                        print txt_5UTR.strip()
                    except:
                        N_5UTR += 1
                        pass

                    # STARTING OF A TRANSCRIPT
                    START = int(temp[3])
                    E_start = 1
                    E_end = 0
                    Length_Trans = 1

            EXON = "exon_" + str(N_Ex)
            trans = TID
            E_start = E_end + 1
            E_end = E_start + (int(temp[4]) - int(temp[3]))

            ########################################
            # Printing CDS Part
            ########################################
            try:
                ###########################################################
                # Gene orientation
                ###########################################################
                if GOrin[trans] == ".":
                    Gene_ori = GOrin_TD[trans]
                else:
                    Gene_ori = GOrin[trans]

                ####################
                # Find 5 UTR Range
                ####################
                #start,end = separate_range(TRANS[TID]["five_prime_UTR"])
                if TRANS[TID].has_key("five_prime_UTR"):
                    start,end = separate_range(TRANS[TID]["five_prime_UTR"])
                    UTR5_start,UTR5_end,txt_5UTR = calculate_5_UTR_regions(start,end,EXONS[trans],gstart,gend,GOrin,GOrin_TD,trans,GeneNames_New,temp)
                else:
                    Tstart,Tend = separate_range(Trans_Range[TID])
                    UTR5_start,UTR5_end = calculate_5_UTR_regions_in_case_of_missing(TRANS[TID],GOrin_TD[trans],Gene_ori,Tstart,Tend)

                #print "5UTR Range:",GeneNames_New[trans],"\t",UTR5_start,"\t",UTR5_end

                #########################################################
                # Irrespective of the strand orientation (+ or -)       #
                # start and end of features are mentioned/considered    #
                # as + oriented to avoid confusion                      #
                #########################################################

                #print GeneNames_New[trans],CDS_START,CDS_END
                #print GeneNames_New[trans],CDS_START,CDS_END
                if TRANS[TID].has_key("three_prime_UTR"):
                    start,end = separate_range(TRANS[TID]["three_prime_UTR"])
                    UTR3_start,UTR3_end,txt_5UTR = calculate_3_UTR_regions(start,end,EXONS[trans],gstart,gend,GOrin,GOrin_TD,trans,GeneNames_New,temp)
                else:
                    #print "3UTR Range Before:",GeneNames_New[trans]
                    Tstart,Tend = separate_range(Trans_Range[TID])
                    #print "3UTR TRange :",GeneNames_New[trans],"\t",Tstart,"\t",Tend
                    UTR3_start,UTR3_end = calculate_3_UTR_regions_in_case_of_missing(TRANS[TID],GOrin_TD[trans],Gene_ori,Tstart,Tend)

                #print "3UTR Range:",GeneNames_New[trans],"\t",UTR3_start,"\t",UTR3_end

                if TRANS[TID].has_key("five_prime_UTR"):
                    UTR5 = 1
                else:
                    UTR5 = 0

                if TRANS[TID].has_key("three_prime_UTR"):
                    UTR3 = 1
                else:
                    UTR3 = 0

                # Find CDS whole region between UTR regions
                if Gene_ori == "+" or Gene_ori == ".":
                    UTR_start = UTR5_end
                    UTR_end = UTR3_start
                else:
                    UTR_start = UTR3_end
                    UTR_end = UTR5_start

                # Avoid confusion
                # Change these silly codes SOON
                CDS_START = UTR_start
                CDS_END = UTR_end

                # Then decide the CDS range
                cds_start,cds_end = separate_range(TRANS[TID]["CDS"])

                mRNA_start,mRNA_end = separate_range(TRANS[TID]["mRNA"])

                CDS = 0
                if GOrin[trans] == ".":
                    GO = GOrin_TD[trans]
                else:
                    GO = GOrin[trans]

                #print "UTR Range:",GeneNames_New[trans],"\t",CDS_START,"\t",CDS_END
                # Check whether the exons fall within the CDS boundary
                #cds_start,cds_end,CDS = decide_on_CDS(cds_start,cds_end,temp[3],temp[4])
                cds_start,cds_end,CDS = decide_on_CDS(CDS_START,CDS_END,temp[3],temp[4],UTR5,UTR3,GO,GeneNames_New[trans])

                #print "CDS Range:",GeneNames_New[trans],"\t",cds_start,"\t",cds_end

                if GOrin[trans] == ".":
                    GO = GOrin_TD[trans]
                else:
                    GO = GOrin[trans]

                if CDS == 1:
                    txt_CDS += scaf + "\t" + "Salmon" + "\t" + "CDS" + "\t" + str(cds_start) + "\t" + str(cds_end) + "\t" + str(temp[5]) + "\t" + GO + "\t" + str(temp[7]) + "\t" + "Parent=" + GeneNames_New[trans] + "\n"
            except:
                pass


# Printing Gene parts
# We are printing because each gene should be printed only once
for i in Gene_Lines:
    print Gene_Lines[i].strip()

# Write the information about the scaffolds where no match transcript was matched
# So that GBrowse can show results to the scaffolds which are mapping in blast result.

