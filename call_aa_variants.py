"""
Created on Fri Aug 23 16:09:03 2019

Call amino acid variants based on predicted effects of SNVs in ref sequences
get_coding_positions defines coding regions and coding frames across TSWV genome

@author: david
"""
from Bio import SeqIO
from Bio.Seq import Seq

def call_aa_variants(variants):
    
    "Get coding positions if calling amino acid variants"
    #if coding_positions is None:
    coding_positions,reverse_sense,aa_positions,annotation = get_coding_positions()
    ref_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/tswv_ref/'
    #ref_fasta = ref_path + 'tswv_ref_full.fasta'
    ref_fasta = ref_path + 'TF2_consensus.fasta'
    seq_dict = SeqIO.to_dict(SeqIO.parse(ref_fasta, "fasta"))
    
    aa_variant = []
    aa_subst = []
    cds_annotation = []
    for index, var in variants.iterrows():
        region = var['REGION']
        pos = var['POS'] - 1
        alt = var['ALT']
        cds = annotation[region][pos]
        ref_seq = seq_dict[region].seq
        nonsyn = False
        if not ("+" in alt or "-" in alt): # ingore indels
            code_pos = coding_positions[region][pos] # position in codon e.g. 1st, 2nd or 3rd position
            reverse = reverse_sense[region][pos] # reverse complement seq?
            if code_pos == 1:
                alt_codon = alt + ref_seq[pos+1:pos+3]
                ref_codon = ref_seq[pos:pos+3]
                nonsyn,alt_aa,ref_aa = isNonSynonomous(alt_codon,ref_codon,reverse)
            elif code_pos == 2:
                alt_codon = Seq(ref_seq[pos-1] + alt + ref_seq[pos+1])
                ref_codon = ref_seq[pos-1:pos+2]
                nonsyn,alt_aa,ref_aa = isNonSynonomous(alt_codon,ref_codon,reverse)
            elif code_pos == 3:
                alt_codon = ref_seq[pos-2:pos] + alt
                ref_codon = ref_seq[pos-2:pos+1]
                nonsyn,alt_aa,ref_aa = isNonSynonomous(alt_codon,ref_codon,reverse)

        subst = ''
        if nonsyn:
            subst = str(ref_aa) + str(aa_positions[region][pos]) + str(alt_aa)
        aa_subst.append(subst)
        aa_variant.append(nonsyn)
        cds_annotation.append(cds)
        
    return aa_variant,aa_subst,cds_annotation

def isNonSynonomous(alt_codon,ref_codon,reverse=False):
    
    nonsyn = False
    if reverse:
        alt_aa = alt_codon.reverse_complement().translate()
        ref_aa = ref_codon.reverse_complement().translate()
    else:
        alt_aa = alt_codon.translate()
        ref_aa = ref_codon.translate()
    if alt_aa != ref_aa:
        nonsyn = True
    return nonsyn,alt_aa,ref_aa
    

def get_coding_positions():
    
    "Get coding infor for seg S"
    # Remember that NSs is encoded as a reverse complement
    segS = [0]*88 # positions 1 - 88 are non-coding
    lengthN = int((1483 - 89 + 1) / 3)
    segS.extend([1,2,3]*lengthN) # positions 89-1483 code for N
    segS.extend([0]*(1986 - 1484 + 1)) # positions 1484 to 1986 are non coding
    lengthNSs = int((2763 - 1987 + 1) / 3) 
    segS.extend([1,2,3]*lengthNSs) # positions 1987-2763 code for NSs
    segS.extend([0]*(2916 - 2764 + 1)) # positions 2764 to 2916 are non coding
    
    rev_segS = [0]*2916
    rev_segS[1987-1:2763] = [1] * (2763-1987+1)
    
    aa_pos_segS = [0]*2916
    aa_pos_segS[89-1:1483] = [x for x in range(1,lengthN+1) for y in range(3)]
    NSs_aa_pos = [x for x in range(1,lengthNSs+1) for y in range(3)]
    NSs_aa_pos.reverse()
    aa_pos_segS[1987-1:2763] = NSs_aa_pos
    
    annotation_segS = ['']*2916
    annotation_segS[89-1:1483] = ['N' for x in range(1,lengthN+1) for y in range(3)]
    annotation_segS[1987-1:2763] = ['NSs' for x in range(1,lengthNSs+1) for y in range(3)]
    
    "Get coding infor for seg M -- treating GP as one protein"
    # Add  NSm 101..1009 and GP 1330..4737
    #segM = [0]*100 # positions 1 - 100 are non-coding
    #lengthNSm = int((1009 - 101 + 1) / 3)
    #segM.extend([1,2,3]*lengthNSm) # positions 101-1009 code for NSm
    #segM.extend([0]*(1329 - 1010 + 1)) # positions 1010 to 1329 are non coding
    #lengthGP = int((4737 - 1330 + 1) / 3)
    #segM.extend([1,2,3]*lengthGP) # positions 1330-4737 code for GP
    #segM.extend([0]*(4821 - 4738 + 1)) # positions 4738 to 4821 are non coding
    
    #rev_segM = [0]*4821
    #rev_segM[1330-1:4737] = [1] * (4737-1330+1)
    
    #aa_pos_segM = [0]*4821
    #aa_pos_segM[101-1:1009] = [x for x in range(1,lengthNSm+1) for y in range(3)]
    #GP_aa_pos = [x for x in range(1,lengthGP+1) for y in range(3)]
    #GP_aa_pos.reverse()
    #aa_pos_segM[1330-1:4737] = GP_aa_pos
    
    #annotation_segM = ['']*4821
    #annotation_segM[101-1:1009] = ['NSm' for x in range(1,lengthNSm+1) for y in range(3)]
    #annotation_segM[1330-1:4737] = ['GP' for x in range(1,lengthGP+1) for y in range(3)]
    
    "Get coding infor for seg M"
    # Add  NSm 101..1009 and GP 1330..4737
    segM = [0]*100 # positions 1 - 100 are non-coding
    lengthNSm = int((1009 - 101 + 1) / 3)
    segM.extend([1,2,3]*lengthNSm) # positions 101-1009 code for NSm
    segM.extend([0]*(1329 - 1010 + 1)) # positions 1010 to 1329 are non coding
    
    lengthGc = int((3285 - 1330 + 1) / 3)
    segM.extend([1,2,3]*lengthGc) # positions 1330-3285 code for Gc   
    lengthGn = int((4737 - 3286 + 1) / 3)    
    segM.extend([1,2,3]*lengthGn) # positions 3286-4737 code for Gn
    segM.extend([0]*(4821 - 4738 + 1)) # positions 4738 to 4821 are non coding
    
    rev_segM = [0]*4821
    rev_segM[1330-1:4737] = [1] * (4737-1330+1)
    
    aa_pos_segM = [0]*4821
    aa_pos_segM[101-1:1009] = [x for x in range(1,lengthNSm+1) for y in range(3)]
    
    Gc_aa_pos = [x for x in range(1,lengthGc+1) for y in range(3)]
    Gc_aa_pos.reverse()
    aa_pos_segM[1330-1:3285] = Gc_aa_pos
    
    Gn_aa_pos = [x for x in range(1,lengthGn+1) for y in range(3)]
    Gn_aa_pos.reverse()
    aa_pos_segM[3286-1:4737] = Gn_aa_pos
    
    annotation_segM = ['']*4821
    annotation_segM[101-1:1009] = ['NSm' for x in range(1,lengthNSm+1) for y in range(3)]
    annotation_segM[1330-1:3285] = ['Gc' for x in range(1,lengthGc+1) for y in range(3)]
    annotation_segM[3286-1:4737] = ['Gn' for x in range(1,lengthGn+1) for y in range(3)]
    
    
    "Get coding infor for seg L"
    # Add RdRp 34..8661
    segL = [0]*33 # positions 1 - 33 are non-coding
    lengthRdRp = int((8661 - 34 + 1) / 3)
    segL.extend([1,2,3]*lengthRdRp) # positions 34-8661 code for RdRp
    segL.extend([0]*(8897 - 8662 + 1)) # positions 8662 to 8897 are non coding
    
    rev_segL = [0]*8897
    
    aa_pos_segL = [0]*8897
    aa_pos_segL[34-1:8661] = [x for x in range(1,lengthRdRp+1) for y in range(3)]
    
    annotation_segL = ['']*8897
    annotation_segL[34-1:8661] = ['RdRp' for x in range(1,lengthRdRp+1) for y in range(3)]
    
    coding_positions = {'TSWV_segS':segS,'TSWV_segM':segM,'TSWV_segL':segL}
    reverse_sense = {'TSWV_segS':rev_segS,'TSWV_segM':rev_segM,'TSWV_segL':rev_segL}
    aa_positions = {'TSWV_segS':aa_pos_segS,'TSWV_segM':aa_pos_segM,'TSWV_segL':aa_pos_segL}
    cds_annotation = {'TSWV_segS':annotation_segS,'TSWV_segM':annotation_segM,'TSWV_segL':annotation_segL}
    
    return coding_positions,reverse_sense,aa_positions,cds_annotation