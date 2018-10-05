import cafa
import constant as const
import taxon

cobj = cafa.cafa()

cobj.create_author_kwd_taxon_dict()
cobj.fmax_kwdscores_by_taxonID(1, 1, 9606)

#cobj.plot_ontology_fmax_results('BPO')
#cobj.plot_ontology_fmax_results('CCO')
cobj.plot_ontology_fmax_results('MFO')


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