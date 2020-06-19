import requests
import json
import os
import re
import gzip
import shutil
import tarfile
from pathlib import Path
import pandas as pd
import numpy as np
from maayanlab_bioinformatics.harmonization import ncbi_genes
import math

def obtain_sample_ids(cancer_type, ensembl = False):
    """
    Obtains sample ids (from cases). In the future, would be better to obtain sample ids 
    while getting the CSVs. If we want ensembl_ids (which requires more labor, including 
    going through each of the downloaded files) then we say ensembl = True. Otherwise,
    it is faster and simpler to just get the sample_ids (doesn't require going through
    text files)
    """
    # Endpoints
    base_url = 'https://api.gdc.cancer.gov/'
    files_endpt = base_url + 'files/'
    genes_endpt = base_url + 'genes/'
    cases_endpt = base_url + 'cases/'
    data_endpt = base_url + "data/"

    # data type of files we want
    data_type = "htseq.counts"

    os.makedirs("TCGA_cancer_downloads", exist_ok = True)
    
    # The 'fields' parameter is passed as a comma-separated string of single names.
    fields = "cases.samples.sample_id,file_id,file_name,cases.case_id"

    # filter files for only RNA-Seq results
    filters = {
        "op": "and",
         "content":[
             {
                "op": "in",
                "content":
                 {
                     "field": "files.experimental_strategy",
                     "value": ["RNA-Seq"],
                 }
             },
             {
                "op": "in",
                "content":
                 {
                     "field": "access",
                     "value": ["open"],

                 }
             },
             {
                "op": "in",
                "content":
                 {
                     "field": "files.file_name",
                     "value": ["*htseq.counts.gz"],
                 }
             },

             {
                "op": "in",
                "content":
                 {
                     "field": "cases.diagnoses.primary_diagnosis",
                     "value": [cancer_type],
                 }
             }
         ],
    }

    # build parameters object
    params = {
        "fields": fields,
        "filters": json.dumps(filters),
        "size": 100000
    }

    # get list of all files with RNA-seq results
    response = requests.get(files_endpt, params = params) # optionally also provide params argument
    data = json.loads(response.content.decode("utf-8"))

    # get list of results
    results = data["data"]["hits"]
    results = list(filter(lambda x: data_type in x["file_name"], results))

    file_uuid_list = [ entry["file_id"] for entry in results]
    case_uuid_list = [ entry["cases"][0]["case_id"] for entry in results]
    sample_id_list = [ entry["cases"][0]["samples"][0]["sample_id"] for entry in results ]
    if ensembl:
        ensembl_id_list = obtain_ensembl_ids(cancer_type, case_uuid_list, file_uuid_list)
        return sample_id_list, ensembl_id_list 
    return sample_id_list

    
def obtain_ensembl_ids(cancer_type, case_uuid_list, file_uuid_list):
    """
    Obtains ensembl ids (the name of genes before they are converted into symbols. 
    In the future, would be easier to obtain these ids while processing.
    """
    
    base_url = 'https://api.gdc.cancer.gov/'
    data_endpt = base_url + "data/"
    
    os.makedirs(f'data/{cancer_type}', exist_ok = True)

    df_files_cases=pd.DataFrame({"case": case_uuid_list },  index=file_uuid_list)
    file_to_case = df_files_cases.to_dict()["case"] # dict mapping file id to case id

    params = {"ids": file_uuid_list}

    # A POST is used, so the filter parameters can be passed directly as a Dict object.
    response = requests.post(data_endpt,data = json.dumps(params),headers={"Content-Type": "application/json"})

    # filename is found in the Content-Disposition header of response
    response_head_cd = response.headers["Content-Disposition"]
    file_name = re.findall("filename=(.+)", response_head_cd)[0]

    downloads_folder = "TCGA_cancer_downloads/"

    # Save .tar.gz zipped file to TCGA_downloads folder
    with open(downloads_folder + file_name, "wb") as f_out:
        f_out.write(response.content)

    # extract the root tar archive
    tar = tarfile.open(downloads_folder + file_name, "r:gz")
    tar.extractall(f'./{downloads_folder}')
    folder = file_name.split(".tar.gz")[0]

    for tarinfo in tar:
        if (tarinfo.name == "MANIFEST.txt"): continue
        file_id = tarinfo.name.split("/")[0]

        # unzip inner .gz files
        with gzip.open(downloads_folder + tarinfo.name, "rb") as f_in:
            os.makedirs(f"data/{cancer_type}",exist_ok = True)
            with open(f"data/{cancer_type}/{file_to_case[file_id]}.txt", "wb") as f_out:
                f_out.write(f_in.read())

    tar.close()

    # initialize empty df
    df = pd.DataFrame({"gene": []})
    df = df.set_index("gene")

    # loop over files, merging with pre-existing data
    for file in Path(f'data/{cancer_type}').glob('*.txt'):
        with open(file, "rb") as f_in:
            new_df = pd.read_csv(f_in, sep = "\t", header = None)
            file_id = re.findall(f"data/{cancer_type}/(.+).txt", f_in.name)[0]
            new_df.columns = ["gene", file_id]
            new_df.gene.replace(to_replace = r'\..*$', value = "", regex=True,
               inplace=True) # collapse all versions of same gene to one gene
            new_df = new_df.set_index("gene")
            df = pd.DataFrame.merge(df, new_df, how="outer", left_on = "gene", right_on = "gene")

    # drop rows not corresponding to genes (i.e. metadata)

    non_genes = list(filter(lambda val: not "ENSG" in val,np.array(df.index.values)))
    df = df.drop(non_genes)

    ncbi = pd.DataFrame(ncbi_genes.ncbi_genes_fetch())

    ensembl_ids = df.index.astype("S").tolist()
    
    return ensembl_ids