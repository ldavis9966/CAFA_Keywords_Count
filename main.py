import cafa
import constant as const
import taxon

cobj = cafa.cafa()


cobj.set_path_cafa_team_files(const.DATA_ROOT_DIRECTORY)
cobj.set_path_zipped_files(const.TAR_DATA_ROOT_DIRECTORY)
cobj.set_path_fmax_files(const.FMAX_FILES_DIRECTORY)

cobj.create_author_kwd_taxon_dict()
cobj.fmax_kwdscores_by_taxonID('BPO', 1, 1, 9606)

print(cobj.keyword_fmax_score)
print(cobj.keyword_relative_fmax_score)

#sum = 0
#for kwd in cobj.keyword_relative_fmax_score:
#    sum += cobj.keyword_relative_fmax_score[kwd]
#print(sum)
#sum = 0
#for kwd in cobj.keyword_fmax_score:
#    sum += cobj.keyword_fmax_score[kwd]
#print(sum)
#print("Taxon:", taxon.taxon_name_converter(7227))