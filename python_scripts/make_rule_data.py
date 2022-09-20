import os
import time
import argparse
import sys
import hashlib
import ujson
import pickle
import yaml

import pandas as pd
import numpy as np
import rubrix as rb

# from io import StringIO
from tqdm import tqdm, trange
from metapub import PubMedFetcher
from Bio import Entrez
sys.path.append('.')
from ids import email, api_key, rubrix_api_key
from yaml import Loader
from collections import defaultdict


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--query_path', help='path to text file with the query to run')
parser.add_argument('--records_logfile', default='logged_records.txt', help='path to list of PMIDs of records that have already been logged')

args = parser.parse_args()


def get_query_hash(x):
    '''
    Create unique hash of a PubMed query
    '''
    return hashlib.sha1(x.encode()).hexdigest()[:50]

def search_pmid(search_string):
    '''
    Example search string: 
    For testing:
    '(("colorectal" AND (neoplasm OR cancer OR tumour)) OR "Colorectal neoplasms"[MeSH]) AND ("Adrenergic beta-antagonists"[MeSH] OR "Antihypertensive Agents"[MeSH] OR "beta-blockers") AND ("Cancer Survivors"[MeSH] OR "cancer survivorship" OR "cancer survivors" OR "cancer survival")'
    
    For usage:
    type_cancer_all = ["colorectal", "lung", "rectum", "female breast", "colon", "ovary", "prostate", "bladder"]
    type_cancer_all_mesh = ["Colorectal", "Lung", "Rectal", "Breast", "Colonic", "Ovarian", "Prostatic", "Urinary bladder"]
    type_cancer = type_cancer[0]
    type_cancer_mesh = type_cancer_all_mesh[0]
    '(("{}" AND (neoplasm OR cancer OR tumour)) OR "{} neoplasms"[MeSH]) AND ("Statins" OR "Hydroxymethylglutaryl-CoA Reductase Inhibitors"[MESH] OR "Adrenergic beta-antagonists"[MeSH] OR "beta-blockers" OR "Antihypertensive Agents"[MeSH] OR "Antihypertensive Agents" OR "Metformin" OR "Metformin"[MESH] OR "Aspirin" OR "Aspirin"[MESH] OR "ARBS" OR "Angiotensin Receptor Antagonists"[MESH] OR "NSAIDs" OR "Anti-Inflammatory Agents, Non-Steroidal"[MESH]) AND ("clinical trails" OR "retrospective")'.format(type_cancer, type_cancer_mesh)
    '''
    fetch = PubMedFetcher(cachedir='./cachedir/')
    pmids = fetch.pmids_for_query(search_string, retmax=5000000)
    return pmids

# def pmid_abstract(pmid):
#     '''
#     Set up email and api_key in utils/ids.py
#     Obtain api key: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
#     '''
#     Entrez.email = email
#     Entrez.api_key = api_key
#     handle = Entrez.efetch(id=pmid, db='pubmed', retmode="xml", rettype='MEDLINE')
#     record = Entrez.read(handle)
#     # print(handle.readline().strip())
#     handle.close()
#     tmp = record['PubmedArticle'][0]['MedlineCitation']['Article']
#     if 'Abstract' in tmp.keys():
#         abstract = ""
#         for part in tmp['Abstract']['AbstractText']:
#             abstract += part
#         return abstract
#     else:
#         return None

def time_msg(fun, inputs, msg):
        t0 = time.time()
        print(msg)
        r = fun(*inputs)
        t1 = time.time()
        print(f"Done: {t1-t0:.2f}s\n")
        return r

def get_journal_info(journal_record):
    record_dict =  {
            # 'issn':journal_record['ISSN'].strip(), 
            'journal_issue':{key:val for key, val in journal_record['JournalIssue'].items()},
            'journal_title':journal_record['Title'],
            'journal_abbreviation': journal_record['ISOAbbreviation']
           }
    if 'ISSN' in journal_record.keys():
        record_dict['issn'] = journal_record['ISSN'].strip()
    

def get_abstracts(pmids):
    '''
    Pull abstracts for a long list of PMIDs

    Inputs:
    --------------------
        pmids: list 
            List of PMIDs to pull from PubMed

    Outputs:
    --------------------
        abstract_dict: dict
            Dictionary of abstracts with appropriate metadata
            Keys of dict:
                'parsed': Abstracts suitable for curation based on exclusion criteria
                'no_substances': Abstracts with text but without chemical substances
                'no_abstract': PMIDs with title only, no abstract
    '''

    # Add in credentials
    Entrez.email = email
    Entrez.api_key = api_key

    # Pull PMIDs in batches
    pmids = set(pmids)
    pmids = [str(x) for x in pmids]
    records = []
    for i in trange(len(pmids)//1000 + 1):
        handle = Entrez.efetch(db='pubmed', id=','.join(pmids[i*1000:(i+1)*1000]), rettype='MEDLINE', retmode='xml')   
        records.extend(Entrez.read(handle)['PubmedArticle'])
        handle.close()
    
    parsed = []
    no_abstract = []
    no_substances = []
    
    for record in tqdm(records):
        # Grab PMID
        pmid = record['MedlineCitation']['PMID'].strip()
        
        # Make sure we have an abstract
        if 'Abstract' in record['MedlineCitation']['Article'].keys():
            abstract = ' '.join(record['MedlineCitation']['Article']['Abstract']['AbstractText'])
            if len(abstract.split()) <= 5:
                no_abstract.append(pmid)
                # print(abstract)
                continue
        else:
            continue
            
        # Get other metadata  
        journal = get_journal_info(record['MedlineCitation']['Article']['Journal'])
        title = record['MedlineCitation']['Article']['ArticleTitle']
        if "AuthorList" in record['MedlineCitation']['Article']:
            authorlist = record['MedlineCitation']['Article']['AuthorList']
        else:
            # print(record['MedlineCitation']['Article'].keys())
            authorlist = []
        
        authors = [x['LastName'] + ', ' + x['ForeName'] for x in authorlist if 'LastName' in x.keys() and 'ForeName' in x.keys()]
        text = title + '\n' + abstract
        
        # Make sure article has some substances
        if 'ChemicalList' in record['MedlineCitation'].keys():
            substance_ui = [x['NameOfSubstance'].attributes['UI'] for x in record['MedlineCitation']['ChemicalList']]
            substances = [x['NameOfSubstance'].strip() for x in record['MedlineCitation']['ChemicalList']]
        else:
            no_substances.append({'pmid':pmid, 'text':text})
            continue
            
        metadata = {'pmid':pmid, 'title': title, 'authors':authors, 'substances':substances, 'substance_mesh_id':substance_ui, 'journal':journal}
        parsed.append({'metadata':metadata, 'abstract':abstract, 'text':text})
        
    return {'parsed':parsed, 'no_abstract':no_abstract, 'no_substances':no_substances}



def load_abstracts(parsed_abstracts, 
                    api_key=rubrix_api_key, 
                    records_logfile=None, 
                    qc_file=None,
                    workspace='spring_2022_rules_stage_1',
                    name=None,
                    min_curators_per_abstract=2,
                    ):
    rb.init("http://localhost:6900/", api_key=api_key)
    rb.set_workspace(workspace)

    # Make sure record hasn't already been assigned
    records = []
    new_logged_records = []
    logged = defaultdict(int)
    if records_logfile is None:
        logged = defaultdict(int)
    elif os.path.isfile(records_logfile):
        logged_records = set(open(records_logfile, 'r').read().strip().split('\n'))
        for record in logged_records:
            logged[record] += 1
    else:
        logged = defaultdict(int)
    
    for record in tqdm(parsed_abstracts):
        # See how many people have been assigned this abstract
        if logged is not None:
            if logged[record['metadata']['pmid']] >= min_curators_per_abstract:
                # if record not in qc_records:
                    continue

        # Appending to the record list
        records.append(rb.TextClassificationRecord(
                inputs=record["text"],
                prediction=[('population_size',0.2), 
                            ('quantitative_effect_measure',0.2),
                            ('study_drug',0.2),
                            ('control_group',0.2),
                            ('target_disease',0.2)],
                metadata=record['metadata'],
                multi_label=False,  # we also need to set the multi_label option in Rubrix
            )
        )
        new_logged_records.append(record['metadata']['pmid'])
    print("Num records: ", len(records))
    # print('\n'.join(new_logged_records))

    if name is None:
        name = workspace
    rb.log(
        records=records,
        name=name,
        # tags={
        #     "task": "text-classification",
        #     "family": "text-classification",
        #     "dataset": "cancer_stage_1",
        # },
    )

    # if records_logfile is not None:
    #     with open(records_logfile, 'a') as f:
    #         for r in new_logged_records:
    #             f.write(r + '\n')



def main():
    all_pmids = []
    
    # Pull list of PMIDs for many cancer types
    type_cancer_all = ["colorectal", "lung", "rectum", "female breast", "colon", "ovary", "prostate", "bladder"]
    type_cancer_all_mesh = ["Colorectal", "Lung", "Rectal", "Breast", "Colonic", "Ovarian", "Prostatic", "Urinary bladder"]
    for cancer_type, mesh_type in zip(type_cancer_all, type_cancer_all_mesh):
        query = '''(("{}" AND (neoplasm OR cancer OR tumour)) OR "{}"[MeSH]) AND ("Statins" OR "Hydroxymethylglutaryl-CoA Reductase Inhibitors"[MESH] OR "Adrenergic beta-antagonists"[MeSH] OR "beta-blockers" OR "Antihypertensive Agents"[MeSH] OR "Antihypertensive Agents" OR "Metformin" OR "Metformin"[MESH] OR "Aspirin" OR "Aspirin"[MESH] OR "ARBS" OR "Angiotensin Receptor Antagonists"[MESH] OR "NSAIDs" OR "Anti-Inflammatory Agents, Non-Steroidal"[MESH]) AND ("clinical trial" OR "retrospective" OR "prospective" OR "case control" OR "case-control")'''.format(cancer_type, mesh_type)

        all_pmids.extend(search_pmid(query))


    # Grab abstracts for all pmids and dump them to files
    abstract_dict = get_abstracts(all_pmids)
    load_abstracts(abstract_dict['parsed'])
    
    # Upload assignments to individual workspaces
    users_list = yaml.load(open('.users.yaml','r'), Loader)
    users = [u['username'] for u in users_list if 'credits' in u]
    name ='rules_test_2'
    load_abstracts(abstract_dict['parsed'],  name=name)
    # for user in users:
    #     workspace = 'rules_stage_1_' + user
    
    #     load_abstracts(abstract_dict['parsed'], workspace=workspace)

if __name__=='__main__':
    main()