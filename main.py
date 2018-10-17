import cafa
import constant as const
import taxon

cobj = cafa.Cafa()

cobj.create_author_kwd_taxon_dict()
#cobj.fmax_kwdscores_by_taxonID(1, 1, 9606)

print(cobj.author_list)

#cobj.fmax_kwdscores_by_taxon('BPO', taxonName='all')
#cobj.plot_ontology_fmax_results('BPO', taxonName='all')
#cobj.csv_ontology_fmax_results('BPO', 1, 1, 'all')

#cobj.fmax_kwdscores_by_taxon('CCO', taxonName='all')
#cobj.plot_ontology_fmax_results('CCO', taxonName='all')
#cobj.csv_ontology_fmax_results('CCO', 1, 1, 'all')

cobj.fmax_kwdscores_by_taxon('MFO', taxonName='all')
cobj.plot_ontology_fmax_results('MFO', taxonName='all')
cobj.csv_ontology_fmax_results('MFO', 1, 1, 'all')


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
cobj.create_taxonID_dictionary()
print(cobj.taxonID_counts)
cobj.create_taxonID_list()
print(cobj.taxonID_list)

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