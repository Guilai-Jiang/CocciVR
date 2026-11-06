# SPDX-License-Identifier: GPL-3.0-or-later

enterotoxin_genes = {'sea' : '1','selq' : '1','selk' : '1','seb' : '1','seh' : '1','sed' : '1',
                     'see' : '1','sell' : '1','sec' : '1','sec3' : '1','yent2' : '1','yent1' : '1',
                     'selu' : '1','selp' : '1','selo' : '1','seln' : '1','selm' : '1','sei' : '1',
                     'seg' : '1','sej' : '1','selr' : '1','selu2' : '1','selv' : '1'}
                    

Blood_infection_genes = {'eta' : '1','etb' : '1',
                         'tsst' : '1','tsst-1' : '1',
                         'scn' : '2',
                         'lukF-PV' : '3','lukS-PV' : '3'}


'''Anti MRSA antibiotics'''
# STREPTOGRAMIN_genes = ['cfr(B)', 'cfr(C)', 'cfr(D)', 'cfr(E)', 'cipA', 'clbA', 'clbB', 'clbC', 'cplR', 
#         'erm(30)', 'erm(31)', 'erm(32)', 'erm(43)', 'erm(44)', 'erm(45)', 'erm(46)', 'erm(47)', 
#         'erm(48)', 'erm(51)', 'erm(56)', 'erm(B)', 'erm(C)', 'erm(D)', 'erm(E)', 'erm(F)', 
#         'erm(G)', 'erm(H)', 'erm(K)', 'erm(N)', 'erm(O)', 'erm(Q)', 'erm(S)', 'erm(T)', 
#         'erm(U)', 'erm(V)', 'erm(W)', 'erm(X)', 'erm(Y)', 'erm(Z)', 'lsa(A)', 'lsa(B)', 
#         'lsa(C)', 'lsa(D)', 'lsa(E)', 'msr(A)', 'msr(C)', 'msr(D)', 'msr(E)', 'msrF', 'msr(G)', 
#         'msrH', 'sal(A)', 'sal(B)', 'sal(C)', 'sal(D)', 'vat(A)', 'vat(B)', 'vat(C)', 'vat(D)', 
#         'vat(E)', 'vat(F)', 'vat(H)', 'vatI', 'vgb(A)', 'vgb(B)', 'vgbC', 'vmlR']

VANCOMYCIN_genes = [
    'murG', 'rpoC', 'vanA', 'vanA-Ao2', 'vanA-Pa', 'vanA-Pt', 'vanA-Pt2', 'vanA-Sc', 'vanB', 'vanC', 
    'vanC1', 'vanD', 'vanE', 'vanF', 'vanG', 'vanH', 'vanH-A', 'vanH-Ac1', 'vanH-Ao1', 'vanH-Ao2', 
    'vanH-B', 'vanH-D', 'vanH-F', 'vanH-M', 'vanH-O', 'vanH-P', 'vanH-Pt', 'vanH-Pt2', 'vanH-Sc', 
    'vanI', 'vanJ', 'vanK-I', 'vanK-Sc', 'vanL', 'vanM', 'vanN', 'vanO', 'vanP', 'vanR-A', 'vanR-B', 
    'vanR-C', 'vanR-Cd', 'vanR-D', 'vanR-E', 'vanR-F', 'vanR-G', 'vanR-I', 'vanR-L', 'vanR-M', 'vanR-N', 
    'vanR-O', 'vanR-P', 'vanR-Sc', 'vanS-A', 'vanS-B', 'vanS-C', 'vanS-Cd', 'vanS-D', 'vanS-E', 
    'vanS-F', 'vanS-G', 'vanS-I', 'vanS-L', 'vanS-M', 'vanS-N', 'vanS-O', 'vanS-P', 'vanS-Pt', 
    'vanS-Pt2', 'vanS-Sc', 'vanTc', 'vanT-C', 'vanT-Cd', 'vanT-E', 'vanT-G', 'vanTm', 'vanT-N', 
    'vanTr', 'vanU-G', 'vanW-B', 'vanW-G', 'vanW-I', 'vanW-Pt', 'vanX-A', 'vanX-Ac1', 'vanX-Ao1', 
    'vanX-Ao2', 'vanX-B', 'vanX-D', 'vanX-F', 'vanX-I', 'vanX-M', 'vanX-O', 'vanX-P', 'vanX-Pt', 
    'vanX-Sc', 'vanXY', 'vanXY-C', 'vanXY-E', 'vanXY-L', 'vanXY-N', 'vanY-A', 'vanY-B', 'vanY-D', 
    'vanY-F', 'vanY-G', 'vanY-G1', 'vanY-M', 'vanY-N', 'vanY-Pt', 'vanZ1', 'vanZ-A', 'vanZ-F', 'vanZ-Pt']

DAPTOMYCIN_genes = ['cdsA', 'cls', 'dltC', 'gshF', 'liaF', 'liaR', 'liaS', 'mprF', 'pgsA', 'rpoB', 'rpoC', 'walK', 'yybT']

LINEZOLID_genes = ['poxtA', 'poxtA-Ef', 'rplC', 'rplD']

TIGECYCLINE_genes = ['acrB', 'acrR', 'adeL', 'adeS', 'lon', 'mepA', 'oqxR', 'ramR', 'rpsJ', 'tet(X2)', 'tet(X3)', 'tet(X4)', 'tet(X5)', 'tmexC', 'tmexD', 'toprJ']

def calculate_virulence(dict_list):
    enterotoxin_score = 0
    Blood_infection_score = 0 
    for dict in dict_list:
        gene = dict['gene']
        if gene in enterotoxin_genes:
            enterotoxin_score += int(enterotoxin_genes[gene])
        if gene in Blood_infection_genes:
            Blood_infection_score += int(Blood_infection_genes[gene])
    #print(enterotoxin_score,Blood_infection_score)
    return [enterotoxin_score,Blood_infection_score]

def calculate_amr(dict_list):
    amr_score = 0
    mec_gene_list = []
    special_gene_set = set()


    for gene_dict in dict_list:
        gene = gene_dict.get('gene', '')  
        if gene in ('mecA', 'mecC'):
            mec_gene_list.append(gene)
            # if gene in STREPTOGRAMIN_genes:
            #         special_gene_set.add('STREPTOGRAMIN')
        if gene in VANCOMYCIN_genes:
            special_gene_set.add('VANCOMYCIN')
        if gene in DAPTOMYCIN_genes:
            special_gene_set.add('DAPTOMYCIN')
        if gene in LINEZOLID_genes:
            special_gene_set.add('LINEZOLID')
        if gene in TIGECYCLINE_genes:
            special_gene_set.add('TIGECYCLINE') 
    
    #print(special_gene_set)

    amr_score = 2 if mec_gene_list and special_gene_set else (1 if mec_gene_list else 0)



    return mec_gene_list,amr_score,special_gene_set