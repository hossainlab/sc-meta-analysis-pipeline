"""
File handling utilities
"""

import os
import h5py
import pandas as pd
import scanpy as sc
from pathlib import Path
from typing import List, Dict, Optional, Union
import logging

def discover_data_files(data_dir: Union[str, Path]) -> Dict[str, List[Path]]:
    """Discover data files in a directory"""
    
    data_path = Path(data_dir)
    
    if not data_path.exists():
        return {}
    
    file_types = {
        'h5ad': list(data_path.glob('*.h5ad')),
        'h5': list(data_path.glob('*.h5')),
        'csv': list(data_path.glob('*.csv')),
        'tsv': list(data_path.glob('*.tsv')),
        'txt': list(data_path.glob('*.txt')),
        'xlsx': list(data_path.glob('*.xlsx')),
        'mtx': list(data_path.glob('*.mtx'))
    }
    
    # Filter empty lists
    file_types = {k: v for k, v in file_types.items() if v}
    
    return file_types

def load_data_file(file_path: Union[str, Path], 
                  file_type: Optional[str] = None) -> sc.AnnData:
    """Load a single data file"""
    
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    # Infer file type if not provided
    if file_type is None:
        file_type = file_path.suffix.lower().lstrip('.')
    
    logger = logging.getLogger(__name__)
    logger.info(f"Loading {file_type} file: {file_path}")
    
    if file_type == 'h5ad':
        adata = sc.read_h5ad(file_path)
    elif file_type == 'h5':
        adata = sc.read_10x_h5(file_path)
    elif file_type in ['csv', 'tsv', 'txt']:
        # Assume genes are rows, cells are columns
        delimiter = '\t' if file_type in ['tsv', 'txt'] else ','
        df = pd.read_csv(file_path, delimiter=delimiter, index_col=0)
        adata = sc.AnnData(X=df.T)  # Transpose so cells are observations
    elif file_type == 'xlsx':
        df = pd.read_excel(file_path, index_col=0)
        adata = sc.AnnData(X=df.T)
    elif file_type == 'mtx':
        # Assume 10X format
        adata = sc.read_10x_mtx(file_path.parent)
    else:
        raise ValueError(f"Unsupported file type: {file_type}")
    
    # Make variable names unique
    adata.var_names_unique()
    
    return adata

def save_results(adata: sc.AnnData, 
                output_dir: Union[str, Path],
                formats: List[str] = None) -> Dict[str, Path]:
    """Save results in multiple formats"""
    
    if formats is None:
        formats = ['h5ad']
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    saved_files = {}
    
    for fmt in formats:
        if fmt == 'h5ad':
            file_path = output_path / 'processed_data.h5ad'
            adata.write_h5ad(file_path)
            saved_files['h5ad'] = file_path
            
        elif fmt == 'csv':
            # Save expression matrix
            expr_path = output_path / 'expression_matrix.csv'
            pd.DataFrame(
                adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
                index=adata.obs_names,
                columns=adata.var_names
            ).to_csv(expr_path)
            
            # Save metadata
            obs_path = output_path / 'cell_metadata.csv'
            adata.obs.to_csv(obs_path)
            
            var_path = output_path / 'gene_metadata.csv'
            adata.var.to_csv(var_path)
            
            saved_files['csv'] = {
                'expression': expr_path,
                'obs': obs_path,
                'var': var_path
            }
            
        elif fmt == 'xlsx':
            with pd.ExcelWriter(output_path / 'results.xlsx') as writer:
                # Expression data (first 1000 most variable genes)
                if adata.n_vars > 1000:
                    top_genes = adata.var.sort_values('dispersions_norm', ascending=False).index[:1000]
                    expr_data = pd.DataFrame(
                        adata[:, top_genes].X.toarray() if hasattr(adata.X, 'toarray') else adata[:, top_genes].X,
                        index=adata.obs_names,
                        columns=top_genes
                    )
                else:
                    expr_data = pd.DataFrame(
                        adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
                        index=adata.obs_names,
                        columns=adata.var_names
                    )
                
                expr_data.to_excel(writer, sheet_name='Expression')
                adata.obs.to_excel(writer, sheet_name='Cell_Metadata')
                adata.var.to_excel(writer, sheet_name='Gene_Metadata')
            
            saved_files['xlsx'] = output_path / 'results.xlsx'
    
    return saved_files

def check_data_integrity(adata: sc.AnnData) -> Dict[str, bool]:
    """Check data integrity"""
    
    checks = {}
    
    # Check for NaN values
    checks['no_nan_in_X'] = not pd.isna(adata.X).any()
    checks['no_nan_in_obs'] = not adata.obs.isna().any().any()
    checks['no_nan_in_var'] = not adata.var.isna().any().any()
    
    # Check for negative values (assuming count data)
    checks['no_negative_values'] = (adata.X >= 0).all()
    
    # Check dimensions
    checks['consistent_dimensions'] = (
        adata.X.shape[0] == len(adata.obs) and 
        adata.X.shape[1] == len(adata.var)
    )
    
    # Check for empty cells/genes
    checks['no_empty_cells'] = (adata.X.sum(axis=1) > 0).all()
    checks['no_empty_genes'] = (adata.X.sum(axis=0) > 0).all()
    
    return checks

def backup_data(adata: sc.AnnData, backup_dir: Union[str, Path]) -> Path:
    """Create a backup of the data"""
    
    backup_path = Path(backup_dir)
    backup_path.mkdir(parents=True, exist_ok=True)
    
    backup_file = backup_path / 'data_backup.h5ad'
    adata.write_h5ad(backup_file)
    
    return backup_file
