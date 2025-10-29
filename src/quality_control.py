"""
Quality Control Module
Implements comprehensive quality control for single-cell RNA-seq data
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from pathlib import Path

class QualityController:
    """Handles quality control filtering for single-cell data"""
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
    def filter_cells_and_genes(self, adata):
        """Main QC filtering function"""
        self.logger.info(f"Starting QC for dataset with {adata.shape[0]} cells, {adata.shape[1]} genes")
        
        # Make a copy to avoid modifying original
        adata = adata.copy()
        
        # Calculate QC metrics
        adata = self._calculate_qc_metrics(adata)
        
        # Plot QC metrics before filtering
        self._plot_qc_metrics(adata, "before_filtering")
        
        # Filter cells
        adata = self._filter_cells(adata)
        
        # Filter genes
        adata = self._filter_genes(adata)
        
        # Plot QC metrics after filtering
        self._plot_qc_metrics(adata, "after_filtering")
        
        self.logger.info(f"After QC: {adata.shape[0]} cells, {adata.shape[1]} genes")
        
        return adata
    
    def _calculate_qc_metrics(self, adata):
        """Calculate quality control metrics"""
        
        # Mitochondrial genes
        adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
        
        # Ribosomal genes
        adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL', 'Rps', 'Rpl'))
        
        # Hemoglobin genes
        adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]', case=False, regex=True)
        
        # Calculate QC metrics
        sc.pp.calculate_qc_metrics(
            adata, 
            percent_top=None, 
            log1p=False, 
            inplace=True,
            var_type=['mt', 'ribo', 'hb']
        )
        
        return adata
    
    def _filter_cells(self, adata):
        """Filter cells based on QC metrics"""
        
        # Get filtering thresholds from config
        min_genes = self.config.get('min_genes_per_cell', 200)
        max_genes = self.config.get('max_genes_per_cell', 5000)
        max_mt_percent = self.config.get('max_mt_percent', 20)
        min_counts = self.config.get('min_counts_per_cell', 1000)
        max_counts = self.config.get('max_counts_per_cell', 30000)
        
        # Filter based on number of genes
        sc.pp.filter_cells(adata, min_genes=min_genes)
        
        # Filter based on total counts
        cell_filter = (
            (adata.obs['total_counts'] >= min_counts) &
            (adata.obs['total_counts'] <= max_counts) &
            (adata.obs['n_genes_by_counts'] <= max_genes) &
            (adata.obs['pct_counts_mt'] <= max_mt_percent)
        )
        
        if 'pct_counts_ribo' in adata.obs.columns:
            max_ribo_percent = self.config.get('max_ribo_percent', 50)
            cell_filter = cell_filter & (adata.obs['pct_counts_ribo'] <= max_ribo_percent)
        
        adata = adata[cell_filter, :].copy()
        
        self.logger.info(f"Filtered cells: {cell_filter.sum()} cells remaining")
        
        return adata
    
    def _filter_genes(self, adata):
        """Filter genes based on expression"""
        
        min_cells = self.config.get('min_cells_per_gene', 3)
        
        # Filter genes expressed in minimum number of cells
        sc.pp.filter_genes(adata, min_cells=min_cells)
        
        self.logger.info(f"Filtered genes: {adata.shape[1]} genes remaining")
        
        return adata
    
    def _plot_qc_metrics(self, adata, suffix=""):
        """Create QC plots"""
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle(f'Quality Control Metrics {suffix}', fontsize=16)
        
        # Total counts per cell
        axes[0, 0].hist(adata.obs['total_counts'], bins=50, alpha=0.7)
        axes[0, 0].set_xlabel('Total counts')
        axes[0, 0].set_ylabel('Number of cells')
        axes[0, 0].set_yscale('log')
        
        # Number of genes per cell
        axes[0, 1].hist(adata.obs['n_genes_by_counts'], bins=50, alpha=0.7)
        axes[0, 1].set_xlabel('Number of genes')
        axes[0, 1].set_ylabel('Number of cells')
        
        # Mitochondrial gene percentage
        axes[0, 2].hist(adata.obs['pct_counts_mt'], bins=50, alpha=0.7)
        axes[0, 2].set_xlabel('Mitochondrial gene %')
        axes[0, 2].set_ylabel('Number of cells')
        
        # Scatter plots
        axes[1, 0].scatter(adata.obs['total_counts'], adata.obs['n_genes_by_counts'], alpha=0.5, s=1)
        axes[1, 0].set_xlabel('Total counts')
        axes[1, 0].set_ylabel('Number of genes')
        
        axes[1, 1].scatter(adata.obs['total_counts'], adata.obs['pct_counts_mt'], alpha=0.5, s=1)
        axes[1, 1].set_xlabel('Total counts')
        axes[1, 1].set_ylabel('Mitochondrial gene %')
        
        axes[1, 2].scatter(adata.obs['n_genes_by_counts'], adata.obs['pct_counts_mt'], alpha=0.5, s=1)
        axes[1, 2].set_xlabel('Number of genes')
        axes[1, 2].set_ylabel('Mitochondrial gene %')
        
        plt.tight_layout()
        
        # Save plot
        output_dir = Path('results/figures')
        output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_dir / f'qc_metrics_{suffix}.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def detect_doublets(self, adata):
        """Detect doublets using scrublet"""
        try:
            import scrublet as scr
            
            scrub = scr.Scrublet(adata.X)
            doublet_scores, predicted_doublets = scrub.scrub_doublets()
            
            adata.obs['doublet_score'] = doublet_scores
            adata.obs['predicted_doublet'] = predicted_doublets
            
            self.logger.info(f"Detected {predicted_doublets.sum()} potential doublets")
            
            return adata
            
        except ImportError:
            self.logger.warning("Scrublet not available, skipping doublet detection")
            return adata
    
    def generate_qc_report(self, adata_before, adata_after):
        """Generate a comprehensive QC report"""
        
        report = {
            'cells_before': adata_before.shape[0],
            'genes_before': adata_before.shape[1],
            'cells_after': adata_after.shape[0],
            'genes_after': adata_after.shape[1],
            'cells_filtered': adata_before.shape[0] - adata_after.shape[0],
            'genes_filtered': adata_before.shape[1] - adata_after.shape[1],
            'filtering_rate_cells': (adata_before.shape[0] - adata_after.shape[0]) / adata_before.shape[0] * 100,
            'filtering_rate_genes': (adata_before.shape[1] - adata_after.shape[1]) / adata_before.shape[1] * 100
        }
        
        # Save report
        output_dir = Path('results/tables')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        report_df = pd.DataFrame([report])
        report_df.to_csv(output_dir / 'qc_report.csv', index=False)
        
        self.logger.info("QC report saved")
        
        return report
