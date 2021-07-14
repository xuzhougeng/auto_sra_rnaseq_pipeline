
import json
import sqlite3

import pandas as pd
import numpy as np

from os.path import basename, join
from urllib import request
from urllib.request import urlopen
from urllib.parse import quote
from urllib.parse import urljoin


pd.set_option("display.max_columns", None)



# calculate the hash value
def myhash(string, size=8):
    import math
    string = string.replace("GSM", "")
    string = string.replace("_","0")[0:100]
    hash_value = int(string)<< 4
    hash_value = int( abs(hash_value) / math.pow(10, size) )
    hash_value = str(abs(hash_value))[0:8]

    return int(hash_value)

# meta_files: metadata sample files
def build_metadata_table(meta_files:list, table_name:str = None, db:str = None):
    
    hash_set = set()
    sample_files = [] 
    hash_values = []
    counts_files = []
    deseq_files = []
    accession = [] 

    # append existed hash and meta
    existed_hash_set = set()
    existed_meta_files = set()
    if table_name is not None and db is not None:
        con = sqlite3.connect(db)
        existed_df = table_from_sql(table_name=table_name, db=db)
        for idx, row in existed_df.iterrows():
            existed_hash_set.add(row['hash'])
            existed_meta_files.add(row['meta_file'])
        con.close()

    for file in meta_files:
        if basename(file) in existed_meta_files:
            continue

        df = pd.read_csv(file, sep = "\t",encoding= 'unicode_escape')
        dict_key = "_".join(sorted(df['GSM'].to_list()))
        hash_value = myhash(dict_key)
        
        # avoid hash collision
        while (hash_value in hash_set or hash_value in existed_hash_set):
            hash_value += 1
        hash_set.add(hash_value)
        hash_values.append(hash_value)
        GSE_ID = np.unique(df['GSE'])[0]
        gene   = np.unique(df['gene'])[0]
        
        file_name = "03_merged_counts/{}_{}_{}.tsv".format(GSE_ID, gene, hash_value)
        deseq_name = "05_DGE_analysis/{}_{}_{}.Rds".format(GSE_ID, gene, hash_value)
        
        sample_files.append(basename(file))
        accession.append( "{}_{}_{}".format(GSE_ID, gene, hash_value) )
        counts_files.append(file_name)
        deseq_files.append(deseq_name)

    df = pd.DataFrame(data = {'hash': hash_values,
                              'accession': accession,
                              'meta_file': sample_files,
                              'count_file': counts_files,
                              'deseq_file': deseq_files,
                              'status': [0 for i in range(len(hash_values)) ]} )
    
    return df

# meta_files: metadata sample files
def build_sample_table(df:pd.DataFrame, file_dir:str = None, sep = "\t"):
    
    sample_dict = {}

    for idx, row in df.iterrows():
        accession = row['accession']
        file_path = join( file_dir, row['meta_file'])
        tmp = pd.read_csv(file_path, sep = sep,encoding= 'unicode_escape')
        sample_dict[accession] = tmp

    samples_df = pd.concat(sample_dict.values(), ignore_index=True)
    rep_len = list(map(len, sample_dict.values()))
    samples_df['accession'] = np.repeat(list(sample_dict.keys()), rep_len )

    return samples_df 

def table_to_sql(df, table_name, db):
    # create or connect to databse
    con = sqlite3.connect(db)
    # save record to database
    df.to_sql(name = table_name, con= con, if_exists="append")
    # close the connection
    con.close()
    return

def update_status( sample_name, table_name, db):
    # create or connect to databse
    con = sqlite3.connect(db)
    # save record to database
    cur = con.cursor()
    cmd = "UPDATE {} SET {} = {} WHERE {} = {}".format(table_name, "status", "1", "meta_file" , basename(sample_name))
    cur.execute(cmd)
    con.commit()
    # close the connection
    con.close()
    return


def table_from_sql( table_name, db):
    con = sqlite3.connect(db)
    df = pd.read_sql('SELECT * FROM {}'.format(table_name), con)
    con.close()
    return df

# notification
def bark_notification(api, contents):
    base_url = api
    content = quote(contents)
    full_url = urljoin(base_url,  content)
    urlopen(full_url)

def feishu_notification(api,contents):
    req =  request.Request(api, method="POST") # this will make the method "POST"
    req.add_header('Content-Type', 'application/json')
    data_dict = {
        "msg_type": "text",
        "content": {"text": "进展报告: " + contents}
    }
    data = json.dumps(data_dict).encode()
    resp = urlopen(req, data = data)
    return  resp

