#!/bin/python

import pandas as pd
import anndata as ad
from pathlib import Path
import argparse
import os
import numpy as np
import gemmapy


def argument_parser():
    parser = argparse.ArgumentParser(description="Preprocess data from GEMMA")
    parser.add_argument("--study_name", type=str, help="Name of the study", default="Velmeshev_et_al.1")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
    
    
def get_sample_meta(samples_raw, samples):
    
    sample_names = [x for x in samples["sample_name"]]
    sample_ids = [x.id for x in samples_raw.data]
    # sample_ids = [x.id for x in samples_raw.data]
    organisms = [x.array_design.taxon.scientific_name.lower().replace(" ", "_") for x in samples_raw.data]
    # need to coerce characterstics to dataframes
    sample_meta = [df for df in samples["sample_characteristics"]]
    # add sample id to each df in sample meta
    sample_meta_updated1 = [df.assign(sample_id=sample_ids[i]) for i, df in enumerate(sample_meta)]
    # add organism to each df in sample meta
    sample_meta_updated2 = [df.assign(organism=organisms[i]) for i, df in enumerate(sample_meta_updated1)]
    # add names back
    sample_meta_updated2 = [df.assign(sample_name=sample_names[i]) for i, df in enumerate(sample_meta_updated2)]
    
    # combine all dfs
    sample_meta_combined = pd.concat(sample_meta_updated2)
    sample_meta_combined.drop_duplicates(subset=["sample_id", "category"], inplace=True)
    sample_meta_df = sample_meta_combined.pivot(index=["sample_id","sample_name","organism"], columns="category",values="value").reset_index()
    organism = list(set(organisms))
    if len(organism) > 1:
        Warning("Warning: Multiple organisms found in the study.")
    return sample_meta_df 


def get_cell_type_assignments(cell_type_assignments):

    combined_df = pd.DataFrame()
    
    for i, assignment in enumerate(cell_type_assignments):
        assignment_name = assignment.name
        assignment_protocol = assignment.protocol.name if assignment.protocol else "NA"
        
        cell_type_ids = assignment.cell_type_ids
        cell_types = assignment.cell_types
        cell_type_value_by_id = {ct.id: ct.value for ct in cell_types}
        cell_type_values = [cell_type_value_by_id[ct_id] if ct_id is not None else 
                            'NA' for ct_id in cell_type_ids]
        
        tempdf = pd.DataFrame({
            assignment_protocol: cell_type_values,
        }) 
        # append new columns to the existing DataFrame, rows are in the same order
        if combined_df.empty:
            combined_df = tempdf
        elif assignment_protocol not in combined_df.columns:
            combined_df = pd.concat([combined_df, tempdf], axis=1)
        # if the column already exists, skip it
        elif assignment_protocol in combined_df.columns:
            print(f"Column {assignment_protocol} already exists in the combined DataFrame. Skipping.")
            continue
    return combined_df 
        

def get_cell_level_characteristics(cell_level_characteristics):
  # get data frame of cell level characteristics for each cell 
  # need to deal with duplicate masks
  combined_df = pd.DataFrame()
  for i, clc in enumerate(cell_level_characteristics):
    cell_type_ids = clc._characteristic_ids
    characteristics = clc.characteristics
    categories = [char._category for char in characteristics]
    if not categories:
        print(f"No categories found for cell level characteristics {i}. Skipping.")
        continue
    characteristic_value_by_id = {ch.id: ch.value for ch in characteristics}
    category_name = categories[0]
    cell_type_values = [characteristic_value_by_id[ct_id] if ct_id is not None else 
                        'NA' for ct_id in cell_type_ids]
    tempdf = pd.DataFrame({
        category_name: cell_type_values,

    })
    if category_name not in combined_df.columns:
        combined_df = pd.concat([combined_df, tempdf], axis=1)
    else:
      # skip if the category already exists
      print(f"Category {category_name} already exists in the combined DataFrame. Skipping.")

  return combined_df 
        
            
def get_cell_level_meta(single_cell_dimension):
    # need a data frame of cell level characteristics, cell type assignments, cell ids, and sample names
    
    cell_ids = single_cell_dimension.data.cell_ids
    sample_names = single_cell_dimension.data.bio_assay_ids
    
    cell_type_assignments = single_cell_dimension.data.cell_type_assignments
    cta_df = get_cell_type_assignments(cell_type_assignments)
    
    # iterate over cell_type_assignments
    cell_level_characteristics = single_cell_dimension.data.cell_level_characteristics
    clc_df = get_cell_level_characteristics(cell_level_characteristics)
    
    # concat row wise
    combined_df = pd.concat([cta_df, clc_df], axis=1)
    combined_df["cell_id"] = cell_ids
    combined_df["sample_id"] = sample_names
    return combined_df
    
    
def main():
    args = argument_parser()
    # use arguments for username and password
    # or use environment variables
    if not os.getenv("GEMMA_USERNAME") or not os.getenv("GEMMA_PASSWORD"):
        raise ValueError("GEMMA_USERNAME and GEMMA_PASSWORD environment variables must be set.")
    username = os.getenv("GEMMA_USERNAME")
    password = os.getenv("GEMMA_PASSWORD")
    client = gemmapy.GemmaPy(auth=[username, password], path='dev')
    study_name = args.study_name
    samples_raw = client.raw.get_dataset_samples(study_name)
    samples = client.get_dataset_samples(study_name, use_processed_quantitation_type=False)
    sample_meta_df = get_sample_meta(samples_raw, samples)
   # sample_meta_df = sample_meta_df.to_csv(f"{study_name}/sample_meta.tsv", sep="\t", index=False)
    

    single_cell_dimension = client.raw.get_dataset_single_cell_dimension(study_name)
    cell_level_meta = get_cell_level_meta(single_cell_dimension)
    
    cell_level_meta = cell_level_meta.merge(sample_meta_df, on=["sample_id"])
    # set index to combination of cell_id and sample_id
    # make string combination of cell_id and sample_id
    cell_level_meta["cell_id"] = cell_level_meta["cell_id"].astype(str)
    cell_level_meta["sample_id"] = cell_level_meta["sample_id"].astype(str)
    cell_level_meta["cell_name"] =  cell_level_meta["sample_id"] + "_" + cell_level_meta["cell_id"]
    cell_level_meta.set_index(["cell_name"], inplace=True)
    # write to gzipped tsv
    cell_level_meta.to_csv(f"{study_name}_cell_level_meta.tsv.gz", sep="\t", compression="gzip")
    
if __name__ == "__main__":
    main()
