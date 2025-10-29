"""
Normalization Module
Implements various normalization strategies for single-cell RNA-seq data
"""

import numpy as np
import scanpy as sc
import logging
from pathlib import Path
import matplotlib.pyplot as plt

class Normalizer:
    """Handles normalization and feature selection for single-cell data"""
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
    def normalize(self, adata):
        """Main normalization function"""
        self.logger.info("Starting normalization...")
        
        # Make a copy to avoid modifying original
        adata = adata.copy()
        
        # Store raw counts
        adata.raw = adata
        
        # Normalize counts
        adata = self._normalize_counts(adata)
        
        # Log transform
        if self.config.get('log_transform', True):
            sc.pp.log1p(adata)
            
        # Find highly variable genes
        adata = self._find_highly_variable_genes(adata)
        
        # Scale data
        if self.config.get('scale', True):
            sc.pp.scale(adata, max_value=10)
            
        self.logger.info("Normalization completed")
        
        return adata
    
    def _normalize_counts(self, adata):
        """Normalize total counts per cell"""
        
        method = self.config.get('method', 'total_count')
        target_sum = self.config.get('target_sum', 1e4)
        
        if method == 'total_count':
            # Total count normalization (CPM-like)
            sc.pp.normalize_total(adata, target_sum=target_sum)
            
        elif method == 'median':
            # Median normalization
            sc.pp.normalize_total(adata, target_sum=None)
            
        elif method == 'scran':
            # scran normalization (requires R/rpy2)
            try:
                import subprocess
                import tempfile
                import os
                
                # This is a simplified version - full implementation would use rpy2
                self.logger.warning("scran normalization not fully implemented, using total count instead")
                sc.pp.normalize_total(adata, target_sum=target_sum)
                
            except Exception as e:
                self.logger.warning(f"scran normalization failed: {e}, using total count instead")
                sc.pp.normalize_total(adata, target_sum=target_sum)
                
        else:
            self.logger.warning(f"Unknown normalization method: {method}, using total count")
            sc.pp.normalize_total(adata, target_sum=target_sum)
            
        return adata
    
    def _find_highly_variable_genes(self, adata):
        """Find highly variable genes for downstream analysis"""
        
        method = self.config.get('hvg_method', 'seurat_v3')
        n_top_genes = self.config.get('n_top_genes', 2000)
        
        if method == 'seurat_v3':
            sc.pp.highly_variable_genes(
                adata,
                flavor='seurat_v3',
                n_top_genes=n_top_genes,
                subset=False
            )
        elif method == 'seurat':
            sc.pp.highly_variable_genes(
                adata,
                flavor='seurat',
                min_mean=0.0125,
                max_mean=3,
                min_disp=0.5,
                subset=False
            )
        elif method == 'cellranger':
            sc.pp.highly_variable_genes(
                adata,
                flavor='cell_ranger',
                n_top_genes=n_top_genes,
                subset=False
            )
        
        # Plot highly variable genes
        self._plot_highly_variable_genes(adata)
        
        # Keep only highly variable genes for downstream analysis
        if self.config.get('subset_hvg', True):
            adata = adata[:, adata.var.highly_variable].copy()
            
        self.logger.info(f"Found {adata.var.highly_variable.sum()} highly variable genes")
        
        return adata
    
    def _plot_highly_variable_genes(self, adata):
        """Plot highly variable genes"""
        
        fig, ax = plt.subplots(figsize=(10, 6))
        sc.pl.highly_variable_genes(adata, ax=ax)
        
        # Save plot
        output_dir = Path('results/figures')
        output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_dir / 'highly_variable_genes.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def calculate_size_factors(self, adata):
        """Calculate size factors for normalization"""
        
        # Simple size factor calculation (total UMI per cell)
        size_factors = np.array(adata.X.sum(axis=1)).flatten()
        size_factors = size_factors / np.median(size_factors)
        
        adata.obs['size_factors'] = size_factors
        
        return adata
    
    def batch_aware_normalization(self, adata, batch_key='batch'):
        """Perform batch-aware normalization"""
        
        if batch_key not in adata.obs.columns:
            self.logger.warning(f"Batch key {batch_key} not found, performing standard normalization")
            return self.normalize(adata)
        
        # Normalize within each batch separately
        normalized_adatas = []
        
        for batch in adata.obs[batch_key].unique():
            batch_adata = adata[adata.obs[batch_key] == batch].copy()
            batch_adata = self.normalize(batch_adata)
            normalized_adatas.append(batch_adata)
        
        # Concatenate back
        adata_normalized = normalized_adatas[0].concatenate(
            normalized_adatas[1:],
            batch_key=batch_key,
            batch_categories=list(adata.obs[batch_key].unique())
        )
        
        return adata_normalized
