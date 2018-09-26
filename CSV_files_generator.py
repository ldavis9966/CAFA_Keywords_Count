import csv
import constant as const
import keyword_counts as kc


# Creates two csv files. First csv file outputs all the authors' taxonIDs and models and shows the keywords used
# The second csv file output the total keywords counts and the total keyword counts by model
def create_keycount_csv(authors_list):
    csvHeader = list(const.METHODOLOGY_KEYWORDS)
    csvHeader.insert(0, "Team")
    csvHeader.insert(1, "Taxon ID")
    csvHeader.insert(2, "Model Number")

    kc.count_kwds(authors_list)

    # Create csv file of keywords used for each team's taxonID and model #
    with open(const.CSV_OUTPUT_DIRECTORY+'/team_model_taxonID_keyword.csv', 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csvHeader)
        writer.writeheader()
        for author in authors_list:
            for taxonID in authors_list[author]:
                for model in authors_list[author][taxonID]:
                    outline = authors_list[author][taxonID][model].copy()
                    outline['Team'] = author
                    outline['Taxon ID'] = taxonID
                    outline['Model Number'] = model
                    writer.writerow(outline)  # end of 1st csv file creation

    # Create csv file for total keyword counts and total keywords count by model #
    csvHeader = list(const.METHODOLOGY_KEYWORDS)
    csvHeader.insert(0, "")
    with open(const.CSV_OUTPUT_DIRECTORY+'/total_keyword_counts.csv', 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csvHeader)
        writer1 = csv.writer(csvfile, delimiter=' ', quotechar="", quoting=csv.QUOTE_NONE)

        writer1.writerow('Total_Keyword_Count')
        writer.writeheader()
        writer.writerow(kc.methodology_keyword_count)

        csvHeader[0] = 'Model Number'
        writer = csv.DictWriter(csvfile, fieldnames=csvHeader)
        writer1.writerow('')
        writer1.writerow('Model_Keyword_Count')
        writer.writeheader()
        for i in ('1', '2', '3'):
            row = kc.model_methodology_keyword_count[i]
            row['Model Number'] = i
            writer.writerow(row)
