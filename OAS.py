import os
import glob
import json
import pandas as pd
from datetime import datetime

os.chdir(".")
extension = 'csv'
all_filenames = [i for i in glob.glob('*.{}'.format(extension))]
print("File processing starts...")

list_of_files = sorted(all_filenames,
                        key = lambda x: os.stat(x).st_size, reverse=True)

print("File sorting ends...")
logf = open("CombineService.log", "w")

Combinedframe = pd.DataFrame()
Counter = 1
for file in list_of_files:
    try:
        sequences = pd.read_csv(file, header=1)
        metadata = ','.join(pd.read_csv(file, nrows=0).columns)
        metadata = json.loads(metadata)
        OAS_info = pd.DataFrame.from_dict(metadata, orient="index").T
        Combinedframe = pd.concat([OAS_info, Combinedframe], axis=0)

        seq_align = pd.DataFrame(sequences['sequence_alignment_aa'])
        seq_align['type'] = "Full-antibody"
        seq_align.columns = ['seq','type']

        fwr_1 = pd.DataFrame(sequences['fwr1_aa'])
        fwr_1['type'] = "FWR-1"
        fwr_1.columns = ['seq','type']

        fwr_2 = pd.DataFrame(sequences['fwr2_aa'])
        fwr_2['type'] = "FWR-2"
        fwr_2.columns = ['seq','type']

        fwr_3 = pd.DataFrame(sequences['fwr3_aa'])
        fwr_3['type'] = "FWR-3"
        fwr_3.columns = ['seq','type']

        fwr_4 = pd.DataFrame(sequences['fwr4_aa'])
        fwr_4['type'] = "FWR-4"
        fwr_4.columns = ['seq','type']

        cdr_1 = pd.DataFrame(sequences['cdr1_aa'])
        cdr_1['type'] = "CDR-1"
        cdr_1.columns = ['seq','type']

        cdr_2 = pd.DataFrame(sequences['cdr2_aa'])
        cdr_2['type'] = "CDR-2"
        cdr_2.columns = ['seq','type']

        cdr_3 = pd.DataFrame(sequences['cdr3_aa'])
        cdr_3['type'] = "CDR-3"
        cdr_3.columns = ['seq','type']

        sequences_filtered = pd.concat([seq_align, fwr_1, fwr_2, fwr_3, fwr_4, cdr_1, cdr_2, cdr_3])
        sequences_filtered.loc[:, 'Species'] = str(OAS_info['Species'][0])
        sequences_filtered.loc[:, 'BType'] = str(OAS_info['BType'][0])
        sequences_filtered['Chain'] = file.split("_")[1]
        sequences_filtered['Isotype'] = file.split("_")[2].split(".")[0]
        sequences_filtered.to_csv(file, header=True)
        logf.write("File {0} finished at: {1}\n".format(file, str(datetime.now())))
        print("File " + str(Counter) + " Finished out of " + str(len(list_of_files)))
        Counter = Counter + 1
    except Exception as e:
        logf.write("Failed at {0}: {1}\n".format(file, str(e)))
    finally:
        pass

Combinedframe.to_csv('Combined_OAS_Stats_cont.csv', index=None, header=True)
