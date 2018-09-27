import cafa
import constant as const
import taxon

cobj = cafa.cafa()


cobj.set_path_cafa_team_files(const.DATA_ROOT_DIRECTORY)
cobj.set_path_zipped_files(const.TAR_DATA_ROOT_DIRECTORY)
cobj.set_path_fmax_files(const.FMAX_FILES_DIRECTORY)

cobj.create_author_kwd_taxon_dict()
cobj.fmax_kwdscores_by_taxonID('BPO', 1, 1, 9606)

print("Taxon:", taxon.taxon_name_converter(7227))