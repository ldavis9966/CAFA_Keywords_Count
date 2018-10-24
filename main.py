import cafa3
import constant as const
import taxon
#import keyword_counts


cobj = cafa3.Cafa3()

#cobj.create_author_kwd_taxon_dict()
#cobj.fmax_kwdscores_by_taxonID(1, 1, 9606)

print(cobj.author_list)

# Test raw keyword counts


# Test raw counts plot and csv
#cobj.create_keycount_csv()

#cobj.count_kwds('1')
#cobj.plot_raw_keyword_counts(8, 6, 8, 8, 12, '5C89C4')

#cobj.count_kwds('2')
#cobj.plot_raw_keyword_counts(8, 6, 8, 8, 12, '5C89C4')

#cobj.count_kwds('2')
#cobj.plot_raw_keyword_counts(4, 3, 6, 6, 8, 'ff00ff')
# End test raw counts plot and csv

# Test relative fmax scores and equal weights scores
#cobj.kwdscores_by_taxon('BPO', taxonName='all')
#cobj.plot_ontology_score_results('BPO', taxonName='all')
#cobj.csv_ontology_score_results('BPO', 1, 1, 'all')

#cobj.kwdscores_by_taxon('MFO', taxonName='all')
#cobj.plot_ontology_score_results('MFO', taxonName='all')
#cobj.csv_ontology_score_results('MFO', 1, 1, 'all')

#cobj.kwdscores_by_taxon('CCO', taxonName='all')
#cobj.plot_ontology_score_results('CCO', taxonName='all')
#cobj.csv_ontology_score_results('CCO', 1, 1, 'all')
# End test relative fmax scores and equal weights scores

# Test relative fmax scores and equal weights scores, specific taxon
#cobj.kwdscores_by_taxon('BPO', taxonID=9606)
#cobj.plot_ontology_score_results('BPO', taxonID=9606)
#cobj.csv_ontology_score_results('BPO', 1, 1, taxonID=9606)

#cobj.kwdscores_by_taxon('MFO', taxonID=9606)
#cobj.plot_ontology_score_results('MFO', taxonID=9606)
#cobj.csv_ontology_score_results('MFO', 1, 1, taxonID=9606)

cobj.kwdscores_by_taxon('CCO', taxonID=9606)
cobj.plot_ontology_score_results('CCO', taxonID=9606)
#cobj.csv_ontology_score_results('CCO', 1, 1, taxonID=9606)
# End test relative fmax scores and equal weights scores, specific taxon


#cobj.fmax_kwdscores_by_taxon('BPO', taxonID=9606)
#cobj.plot_ontology_fmax_results('BPO', taxonID=9606)
#cobj.csv_ontology_fmax_results('BPO', 1, 1, 9606)

#cobj.fmax_kwdscores_by_taxon('CCO', taxonID=9606)
#cobj.plot_ontology_fmax_results('CCO', taxonID=9606)

#cobj.fmax_kwdscores_by_taxon('MFO', taxonID=9606)
#cobj.plot_ontology_fmax_results('MFO', taxonID=9606)

#cobj.plot_ontology_fmax_results('CCO')
#cobj.plot_ontology_fmax_results('MFO')
#cobj.csv_ontology_fmax_results(1, 1, 9606)

print("\n\nAuthor List:")
#print(cobj.keyword_relative_fmax_score)
#print(cobj.keyword_fmax_score_template)
#print(cobj.author_list)
cobj.create_taxon_id_dictionary()
print(cobj.taxon_id_counts)
cobj.create_taxon_id_list()
print(cobj.taxon_id_list)

#ontologies_list = ('BPO', 'CCO', 'MFO')
#bpo_taxons = ("ARATH", "BACSU", "CANAX", "DANRE", "DICDI", "DROME", "ECOLI", "eukarya", "HUMAN", "MOUSE", "prokarya", "RAT", "SCHPO")
#cco_taxons = ("ARATH", "CANAX", "DICDI", "DROME", "ECOLI", "eukarya", "HUMAN", "MOUSE", "prokarya", "RAT")
#mfo_taxons = ("ARATH", "BACSU", "DROME", "ECOLI", "eukarya", "HUMAN", "MOUSE", "prokarya", "RAT", "SALTY")

#cobj.fmax_kwdscores_by_ontology_taxon_list(1, 1, ontologies=ontologies_list, bpo_taxons = bpo_taxons, cco_taxons = cco_taxons, mfo_taxons=mfo_taxons)





#DEBUG CODE
#print('TAXON NAME: '+cobj.taxon_name)
#print('TAXON NAME: '+cobj.taxon_name)

#print(cobj.keyword_fmax_score)
#print(cobj.keyword_relative_fmax_score)

#sum = 0
#for kwd in cobj.keyword_relative_fmax_score:
#    sum += cobj.keyword_relative_fmax_score[kwd]
#print(sum)
#sum = 0
#for kwd in cobj.keyword_fmax_score:
#    sum += cobj.keyword_fmax_score[kwd]
#print(sum)
#print("Taxon:", taxon.taxon_name_converter(7227))
# END DEBUG CODE