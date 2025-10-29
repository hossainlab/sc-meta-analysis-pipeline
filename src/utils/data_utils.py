"""
Data utilities for handling various data formats
"""

import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from pathlib import Path
import logging

def load_h5ad(file_path):
    """Load AnnData object from h5ad file"""
    return ad.read_h5ad(file_path)

def load_csv_as_anndata(file_path, delimiter=','):
    """Load CSV file as AnnData object"""
    
    df = pd.read_csv(file_path, delimiter=delimiter, index_col=0)
    
    # Assume genes are columns, cells are rows
    adata = ad.AnnData(X=df.values, obs=pd.DataFrame(index=df.index), 
                       var=pd.DataFrame(index=df.columns))
    
    # Make variable names unique
    adata.var_names_make_unique()
    
    return adata

def load_mtx_as_anndata(mtx_path, features_path, barcodes_path):
    """Load Matrix Market format as AnnData object"""
    
    from scipy.io import mmread
    
    # Load matrix
    X = mmread(mtx_path).T.tocsr()  # Transpose to genes x cells
    
    # Load features/genes
    features = pd.read_csv(features_path, delimiter='	', header=None)
    gene_names = features[1] if features.shape[1] > 1 else features[0]
    
    # Load barcodes/cells
    barcodes = pd.read_csv(barcodes_path, delimiter='	', header=None)[0]
    
    # Create AnnData object
    adata = ad.AnnData(X=X.T, obs=pd.DataFrame(index=barcodes),
                       var=pd.DataFrame(index=gene_names))
    
    adata.var_names_make_unique()
    
    return adata

def save_anndata(adata, file_path):
    """Save AnnData object to h5ad file"""
    adata.write_h5ad(file_path)

def merge_anndata_objects(adatas, batch_key='batch'):
    """Merge multiple AnnData objects"""
    
    if len(adatas) == 1:
        return list(adatas.values())[0]
    
    # Get list of AnnData objects and their names
    adata_list = list(adatas.values())
    batch_names = list(adatas.keys())
    
    # Concatenate
    adata_merged = adata_list[0].concatenate(
        adata_list[1:],
        batch_key=batch_key,
        batch_categories=batch_names
    )
    
    return adata_merged

def filter_genes_by_expression(adata, min_cells=3):
    """Filter genes by minimum number of expressing cells"""
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata

def filter_cells_by_expression(adata, min_genes=200, max_genes=5000):
    """Filter cells by gene expression criteria"""
    sc.pp.filter_cells(adata, min_genes=min_genes)
    
    if max_genes is not None:
        adata = adata[adata.obs.n_genes_by_counts < max_genes, :].copy()
    
    return adata

def add_qc_metrics(adata):
    """Add quality control metrics to AnnData object"""
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    # Add mitochondrial gene percentage
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, 
                               log1p=False, inplace=True)
    
    return adata

def normalize_gene_names(adata, species='human'):
    """Normalize gene names to standard format"""
    
    if species.lower() == 'human':
        # Convert to uppercase (standard for human genes)
        adata.var_names = adata.var_names.str.upper()
    elif species.lower() == 'mouse':
        # Convert to title case (standard for mouse genes)
        adata.var_names = adata.var_names.str.title()
    
    # Make unique
    adata.var_names_make_unique()
    
    return adata

def subset_to_common_genes(adatas):
    """Subset AnnData objects to common genes"""
    
    # Find common genes
    common_genes = set(list(adatas.values())[0].var_names)
    for adata in list(adatas.values())[1:]:
        common_genes = common_genes.intersection(set(adata.var_names))
    
    common_genes = list(common_genes)
    
    # Subset all datasets
    subset_adatas = {}
    for name, adata in adatas.items():
        subset_adatas[name] = adata[:, common_genes].copy()
    
    return subset_adatas, common_genes

def create_pseudobulk(adata, group_by, aggregation='mean'):
    """Create pseudobulk data by aggregating cells"""
    
    pseudobulk_data = []
    group_names = []
    
    for group in adata.obs[group_by].unique():
        group_cells = adata[adata.obs[group_by] == group]
        
        if aggregation == 'mean':
            group_data = np.array(group_cells.X.mean(axis=0)).flatten()
        elif aggregation == 'sum':
            group_data = np.array(group_cells.X.sum(axis=0)).flatten()
        elif aggregation == 'median':
            group_data = np.median(group_cells.X.toarray() if hasattr(group_cells.X, 'toarray') 
                                  else group_cells.X, axis=0)
        
        pseudobulk_data.append(group_data)
        group_names.append(group)
    
    # Create DataFrame
    pseudobulk_df = pd.DataFrame(pseudobulk_data, 
                                 index=group_names,
                                 columns=adata.var_names)
    
    return pseudobulk_df
