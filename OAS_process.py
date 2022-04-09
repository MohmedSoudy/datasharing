import os
import glob
import json
import pandas as pd
from datetime import datetime

from Bio.SeqUtils import seq1

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
        sequences = pd.read_csv(file)
        Variable_region = sequences[sequences.type == "Full-antibody"]
        Variable_region.drop('type', inplace=True, axis=1)
        Variable_region.to_csv(file, header=True, index=False)
        logf.write("File {0} finished at: {1}\n".format(file, str(datetime.now())))
        print("File " + str(Counter) + " Finished out of " + str(len(list_of_files)))
        Counter = Counter + 1
    except Exception as e:
        logf.write("Failed at {0}: {1}\n".format(file, str(e)))
    finally:
        pass

