
def taxon_name_converter(taxonID):
    #convert from taxonomy ID to name (i.e. from 9606 to HUMAN）
    if isinstance(taxonID, int):
        taxonID = str(taxonID)

    taxonTable = {'10116':'RAT','9606':'HUMAN','3702':'ARATH','7955':'DANRE','44689':'DICDI',
    '7227':'DROME','83333':'ECOLI','10090':'MOUSE','208963':'PSEAE',
    '237561':'CANAX','559292':'YEAST','284812':'SCHPO','8355':'XENLA','224308':'BACSU',
    '99287':'SALTY','243232':'METJA','321314':'SALCH','160488':'PSEPK','223283':'PSESM',
    '85962':'HELPY','243273':'MYCGE','170187':'STRPN','273057':'SULSO','all':'all','prokarya':'prokarya','eukarya':'eukarya'}
    return taxonTable[taxonID]

def taxon_ID_converter(taxon_name):
    #convert from taxonomy ID to name (i.e. from 9606 to HUMAN）

    taxonTable = {'RAT': '10116', 'HUMAN': '9606','ARATH': '3702', 'DANRE': '7955','DICDI': '44689',
    'DROME': '7227', 'ECOLI': '83333', 'MOUSE': '10090', 'PSEAE': '208963',
    'CANAX': '237561','YEAST': '559292', 'SCHPO': '284812', 'XENLA': '8355', 'BACSU': '224308',
    'SALTY': '99287', 'METJA': '243232', 'SALCH': '321314', 'PSEPK': '160488', 'PSESM': '223283',
    'HELPY': '85962', 'MYCGE': '243273', 'STRPN': '170187', 'SULSO': '273057','all':'all','prokarya':'prokarya','eukarya':'eukarya'}
    return taxonTable[taxon_name]
