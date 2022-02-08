import os
import time
import torch
import argparse
import sys
import requests
import hashlib
import ujson
import pickle
import yaml

import pandas as pd
import numpy as np
import rubrix as rb

from io import StringIO
from tqdm import tqdm, trange
from metapub import PubMedFetcher
from Bio import Entrez
from ids import email, api_key
from yaml import Loader
from collections import defaultdict
sys.path.append('.')

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--query_path', help='path to text file with the query to run')
parser.add_argument('--records_logfile', default='logged_records.txt', help='path to list of PMIDs of records that have already been logged')

args = parser.parse_args()


def get_query_hash(x):
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
    # Add in credentials
    Entrez.email = email
    Entrez.api_key = api_key

    # Pull PMIDs in batches
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
# def get_pubtator_abstracts(query, format='pubtator', concepts=['disease','chemical']):
#     '''
#     Get PMIDs for a particular query and store results in a unique location
#     '''
#     pmids = time_msg(search_pmid, [query], 'Running query')
#     # print(pmids)
    
#     abstracts = []
    
#     for i in range(len(pmids)//1000 + 1):
#         start = i*1000
#         end = start + 1000
#         json = {'pmids': pmids[start:end], 'concepts':concepts}
#         r = requests.post("https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/"+format , json=json)
#         abstracts.extend(r.text.split('\n\n'))
        
#     return [x for x in abstracts if x != '']

def distribute_abstracts(users, abstract_list):
    '''
    Distribute abstracts to individual curators based on number of credit hours
    '''
    # Calculate total number of credits 
    total_credits = sum([x['credits'] for x in users if 'credits' in x.keys()])
    n = len(abstract_list)
    abs_per_credit = n//total_credits
    chunks_with_extra = n % total_credits

    # Break abstract list into equally sized chunks
    chunks = []
    for i in range(total_credits):
        chunks.append(abstract_list[i*abs_per_credit : (i+1) * abs_per_credit])
        if i < chunks_with_extra:
            chunks[-1].append(abstract_list[-i])

    # Assign each chunk to 2 users
    # Users end up with 2 chunks per credit hour
    iters = 0
    assignments = {}
    for i, u in enumerate(users):
        if 'credits' not in u.keys():
            continue
        # workspace = 'cancer_stage_1_' + u['username']
        end = iters + 2 * u['credits'] % total_credits
        if end < iters:
            user_chunks = chunks[iters:] + chunks[:end]
        else:
            user_chunks = chunks[iters : iters + 2*u['credits']]
        assignments[u['username']] = [abstract for chunk in user_chunks for abstract in chunk]
        iters += 2 * u['credits'] 
        iters %= total_credits
        
    return assignments
        



def load_abstracts(parsed_abstracts, 
                    api_key='615d255a-45d4-4d0f-bc1c-f998f8e91f40', 
                    records_logfile=None, 
                    qc_file=None,
                    workspace='cancer_stage_1_test',
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

    # Load files needed for quality control
    if qc_file is None:
        qc_records = set([])
    elif os.path.isfile(qc_file):
        qc_records = set(open(records_logfile, 'r').read().strip().split('\n'))
    else:
        qc_records = set([])
    
    for record in tqdm(parsed_abstracts):
        # See how many people have been assigned this abstract
        if logged is not None:
            if logged[record['metadata']['pmid']] >= min_curators_per_abstract:
                if record not in qc_records:
                    continue

        # Appending to the record list
        records.append(rb.TextClassificationRecord(
                inputs=record["text"],
                prediction=[('population_size',0), 
                            ('quantitative_effect_measure',0),
                            ('study_drug',0),
                            ('control_group',0),
                            ('target_disease',0)],
                metadata=record['metadata'],
                multi_label=True,  # we also need to set the multi_label option in Rubrix
            )
        )
        # rb.log(
        #     rb.TextClassificationRecord(
        #         inputs=record["text"],
        #         metadata=record['metadata'],
        #         multi_label=True,  # we also need to set the multi_label option in Rubrix
        #     ),
        #     name="cancer_stage_0",
        #     tags={
        #     "task": "multilabel-text-classification",
        #     "family": "text-classification",
        #     "dataset": "cancer_stage_0",
        #     },
        # )
        new_logged_records.append(record['metadata']['pmid'])

    # print('\n'.join(new_logged_records))
    rb.log(
        records=records,
        name=workspace,
        tags={
            "task": "multilabel-text-classification",
            "family": "text-classification",
            "dataset": "cancer_stage_1",
        },
    )

    with open(records_logfile, 'a') as f:
        for r in new_logged_records:
            f.write(r + '\n')



def main():
    all_pmids = []
    # Pull list of PMIDs

    type_cancer_all = ["colorectal", "lung", "rectum", "female breast", "colon", "ovary", "prostate", "bladder"]
    type_cancer_all_mesh = ["Colorectal", "Lung", "Rectal", "Breast", "Colonic", "Ovarian", "Prostatic", "Urinary bladder"]
    for cancer_type, mesh_type in zip(type_cancer_all, type_cancer_all_mesh):
        query = '''(("{}" AND (neoplasm OR cancer OR tumour)) OR "{}"[MeSH]) AND ("Statins" OR "Hydroxymethylglutaryl-CoA Reductase Inhibitors"[MESH] OR "Adrenergic beta-antagonists"[MeSH] OR "beta-blockers" OR "Antihypertensive Agents"[MeSH] OR "Antihypertensive Agents" OR "Metformin" OR "Metformin"[MESH] OR "Aspirin" OR "Aspirin"[MESH] OR "ARBS" OR "Angiotensin Receptor Antagonists"[MESH] OR "NSAIDs" OR "Anti-Inflammatory Agents, Non-Steroidal"[MESH]) AND ("clinical trial" OR "retrospective" OR "prospective" OR "case control" OR "case-control")'''.format(cancer_type, mesh_type)

        all_pmids.extend(search_pmid(query))

    with open('pmid_list.json', 'w') as f:
        ujson.dump(all_pmids, f)

    abstract_dict = get_abstracts(all_pmids)
    with open('parsed.json', 'w') as f:
        ujson.dump(abstract_dict['parsed'], f)

    with open('no_substances.json', 'w') as f:
        ujson.dump(abstract_dict['no_substances'], f)

    # Distribute assignments to users
    users = yaml.load(open('.users.yaml','r'), Loader)
    assigned = distribute_abstracts(users, abstract_dict['parsed'])
    # for key, val in assigned.items():
    #     print(key, len(val))

    # Upload assignments to individual workspaces
    for user, user_abstracts in assigned.items():
        workspace = 'cancer_stage_1_' + user
        # workspace = 'cancer_stage_1_test'
        load_abstracts(user_abstracts, records_logfile=args.records_logfile, workspace=workspace)

if __name__=='__main__':
    main()