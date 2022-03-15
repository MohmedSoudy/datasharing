from bs4 import BeautifulSoup
import re
import os
import wget
import pandas as pd
from datetime import datetime
from urlextract import URLExtract

with open("vdj-IReceptor download.html", 'r', encoding="utf8") as f:
    contents = f.read()

soup = BeautifulSoup(contents)
Study_divs = soup.find_all("div", {"class": "community-project"})

logf = open("CombineService.log", "w")

Combinedframe = pd.DataFrame()

for i in range(len(Study_divs)):
    Title_header = str(Study_divs[i].find_all("h2", {"class": "community-project-title"}))
    Title = re.search(r'"#">(.*?)</', Title_header).group(1)

    Metadata = Study_divs[i].find_all("div", {"class": "col community-metadata"})
    type(Metadata)

    Study_Metadata = ""
    for data in Metadata:
        Study_Metadata += data.text.replace('\n\n', '').replace('                ', '')

    Study_Metadata = re.sub('Grants:.*?\n\n', '', Study_Metadata, flags=re.DOTALL)
    Study_Metadata = Study_Metadata.replace('\n\n', '')

    Metadata_splitted = Study_Metadata.split('\n')
    Header = []
    Data = []
    for entry in Metadata_splitted:
        # Check if entry has : in it
        Count_colon = [pos for pos, char in enumerate(entry) if char == ":"]
        if len(Count_colon) > 1 and "Study ID" in entry:
            Count_colon.pop(0)
            for j in Count_colon:
                entry = entry[:j] + "_" + entry[j + 1:]
        try:
            Header.append(entry.split(":")[0])
            Data.append(entry.split(":")[1])
        except:
            pass

    MetaData_df = pd.DataFrame(list(zip(Header, Data))).T
    MetaData_df.columns = MetaData_df.iloc[0]
    MetaData_df = MetaData_df.drop(MetaData_df.index[0])
    MetaData_df['Title'] = Title
    # Download the files
    extractor = URLExtract()
    Links = set(extractor.find_urls(str(Study_divs[i])))
    # Create folder for the study
    # mode
    mode = 0o777
    Current_path = os.getcwd()
    Study_id = MetaData_df['Study ID'].iloc[0].strip()
    Out_path = Current_path + '/' + Study_id
    try:
        os.mkdir(Out_path, mode)
        Inx = 1
        for Link in Links:
            if "pubmed" in Link:
                continue
            Out_file = Out_path + "/" + Study_id + "_" + str(Inx) + ".zip"
            #wget.download(Link, out=Out_file)
            Inx += 1
    except:
        pass
    parser = pd.io.parsers.base_parser.ParserBase({'usecols': None})
    Combinedframe = pd.concat([df.set_axis(parser._maybe_dedup_names(df.columns), axis=1) for df in [Combinedframe, MetaData_df]], ignore_index=True)
    logf.write("File {0} finished at: {1}\n".format(Study_id, str(datetime.now())))
    print("File {0} finished at: {1}\n".format(Study_id, str(datetime.now())))

Combinedframe.to_csv('Combined_IReceptor_Stats.csv', index=None, header=True)