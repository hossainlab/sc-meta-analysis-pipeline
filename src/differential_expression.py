"""
Differential Expression Module
Implements various methods for finding differentially expressed genes
"""

import numpy as np
import pandas as pd
import scanpy as sc
import logging
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

class DifferentialExpressionAnalyzer:
    """Handles differential expression analysis"""
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
    def find_markers(self, adata):
        """Main function to find marker genes and perform DE analysis"""
        self.logger.info("Starting differential expression analysis...")
        
        results = {}
        
        # Find cell type markers
        if 'cell_type' in adata.obs.columns:
            results['cell_type_markers'] = self._find_cell_type_markers(adata)
        
        # Condition-based DE analysis
        if 'condition' in adata.obs.columns:
            results['condition_de'] = self._find_condition_de_genes(adata)
        
        # Treatment-based DE analysis
        if 'treatment' in adata.obs.columns:
            results['treatment_de'] = self._find_treatment_de_genes(adata)
        
        # Save results
        self._save_de_results(results)
        
        # Plot results
        self._plot_de_results(adata, results)
        
        self.logger.info("Differential expression analysis completed")
        
        return results
    
    def _find_cell_type_markers(self, adata):
        """Find marker genes for each cell type"""
        
        method = self.config.get('method', 'wilcoxon')
        
        # Use scanpy for fast marker gene finding
        sc.tl.rank_genes_groups(
            adata,
            groupby='cell_type',
            method=method,
            use_raw=True if adata.raw is not None else False,
            n_genes=100
        )
        
        # Extract results
        result = adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        
        markers_dict = {}
        for group in groups:
            markers_dict[group] = pd.DataFrame({
                'gene': result['names'][group],
                'logfoldchanges': result['logfoldchanges'][group],
                'pvals': result['pvals'][group],
                'pvals_adj': result['pvals_adj'][group],
                'scores': result['scores'][group]
            })
        
        return markers_dict
    
    def _find_condition_de_genes(self, adata):
        """Find DE genes between conditions"""
        
        method = self.config.get('method', 'wilcoxon')
        
        # Perform DE analysis for each cell type
        de_results = {}
        
        for cell_type in adata.obs['cell_type'].unique():
            cell_subset = adata[adata.obs['cell_type'] == cell_type].copy()
            
            if 'condition' in cell_subset.obs.columns and len(cell_subset.obs['condition'].unique()) > 1:
                
                sc.tl.rank_genes_groups(
                    cell_subset,
                    groupby='condition',
                    method=method,
                    use_raw=True if cell_subset.raw is not None else False
                )
                
                result = cell_subset.uns['rank_genes_groups']
                groups = result['names'].dtype.names if hasattr(result['names'], 'dtype') else result['names'].keys()
                
                for group in groups:
                    key = f"{cell_type}_{group}_vs_rest"
                    de_results[key] = pd.DataFrame({
                        'gene': result['names'][group],
                        'logfoldchanges': result['logfoldchanges'][group],
                        'pvals': result['pvals'][group],
                        'pvals_adj': result['pvals_adj'][group],
                        'scores': result['scores'][group]
                    })
        
        return de_results
    
    def _find_treatment_de_genes(self, adata):
        """Find DE genes between treatments"""
        
        method = self.config.get('method', 'wilcoxon')
        
        # Similar to condition DE but for treatments
        de_results = {}
        
        for cell_type in adata.obs['cell_type'].unique():
            cell_subset = adata[adata.obs['cell_type'] == cell_type].copy()
            
            if 'treatment' in cell_subset.obs.columns and len(cell_subset.obs['treatment'].unique()) > 1:
                
                sc.tl.rank_genes_groups(
                    cell_subset,
                    groupby='treatment',
                    method=method,
                    use_raw=True if cell_subset.raw is not None else False
                )
                
                result = cell_subset.uns['rank_genes_groups']
                groups = result['names'].dtype.names if hasattr(result['names'], 'dtype') else result['names'].keys()
                
                for group in groups:
                    key = f"{cell_type}_{group}_treatment"
                    de_results[key] = pd.DataFrame({
                        'gene': result['names'][group],
                        'logfoldchanges': result['logfoldchanges'][group],
                        'pvals': result['pvals'][group],
                        'pvals_adj': result['pvals_adj'][group],
                        'scores': result['scores'][group]
                    })
        
        return de_results
    
    def _pseudobulk_deseq2(self, adata, groupby, condition):
        """Perform pseudobulk DESeq2 analysis"""
        
        try:
            # This would require R/rpy2 integration
            # For now, we'll use a simplified approach
            self.logger.warning("DESeq2 pseudobulk analysis not fully implemented")
            return {}
            
        except Exception as e:
            self.logger.error(f"DESeq2 analysis failed: {e}")
            return {}
    
    def _scvi_de_analysis(self, adata, groupby):
        """Perform DE analysis using scvi-tools"""
        
        try:
            import scvi
            
            # Setup scVI model
            scvi.model.SCVI.setup_anndata(adata, batch_key=None)
            model = scvi.model.SCVI(adata)
            model.train()
            
            # Perform differential expression
            de_df = model.differential_expression(groupby=groupby)
            
            return de_df
            
        except Exception as e:
            self.logger.error(f"scVI DE analysis failed: {e}")
            return {}
    
    def _save_de_results(self, results):
        """Save DE results to files"""
        
        output_dir = Path('results/tables')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        for analysis_type, analysis_results in results.items():
            if isinstance(analysis_results, dict):
                # Multiple comparisons
                for comparison, df in analysis_results.items():
                    if isinstance(df, pd.DataFrame):
                        filename = f"de_{analysis_type}_{comparison}.csv"
                        df.to_csv(output_dir / filename, index=False)
            elif isinstance(analysis_results, pd.DataFrame):
                # Single comparison
                filename = f"de_{analysis_type}.csv"
                analysis_results.to_csv(output_dir / filename, index=False)
        
        self.logger.info("DE results saved")
    
    def _plot_de_results(self, adata, results):
        """Plot DE analysis results"""
        
        # Create figure directory
        output_dir = Path('results/figures')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Plot cell type markers if available
        if 'cell_type_markers' in results:
            self._plot_marker_genes(adata, results['cell_type_markers'])
        
        # Plot volcano plots for condition DE
        if 'condition_de' in results:
            self._plot_volcano_plots(results['condition_de'])
        
        # Plot heatmaps
        self._plot_de_heatmaps(adata, results)
    
    def _plot_marker_genes(self, adata, markers_dict):
        """Plot marker genes for cell types"""
        
        # Dotplot of top marker genes
        top_genes = []
        for cell_type, markers in markers_dict.items():
            if isinstance(markers, pd.DataFrame) and len(markers) > 0:
                top_genes.extend(markers.head(5)['gene'].tolist())
        
        if len(top_genes) > 0:
            # Remove duplicates while preserving order
            top_genes = list(dict.fromkeys(top_genes))
            
            # Create dotplot
            fig, ax = plt.subplots(figsize=(12, 8))
            
            try:
                sc.pl.dotplot(
                    adata, 
                    top_genes[:50],  # Limit to top 50 to avoid overcrowding
                    groupby='cell_type',
                    ax=ax,
                    show=False
                )
                plt.title('Top Marker Genes by Cell Type')
                plt.tight_layout()
                plt.savefig(Path('results/figures') / 'marker_genes_dotplot.png', 
                           dpi=300, bbox_inches='tight')
                plt.close()
            except Exception as e:
                self.logger.warning(f"Could not create dotplot: {e}")
                plt.close()
    
    def _plot_volcano_plots(self, de_results):
        """Create volcano plots for DE results"""
        
        for comparison, df in list(de_results.items())[:4]:  # Limit to first 4 comparisons
            if isinstance(df, pd.DataFrame) and len(df) > 0:
                
                fig, ax = plt.subplots(figsize=(10, 8))
                
                # Calculate -log10(p-value)
                df = df.copy()
                df['-log10_pval'] = -np.log10(df['pvals_adj'].clip(lower=1e-300))
                
                # Color points
                colors = np.where(
                    (df['pvals_adj'] < 0.05) & (np.abs(df['logfoldchanges']) > 0.5),
                    'red', 'gray'
                )
                
                ax.scatter(df['logfoldchanges'], df['-log10_pval'], 
                          c=colors, alpha=0.6, s=10)
                
                ax.set_xlabel('Log2 Fold Change')
                ax.set_ylabel('-Log10 Adjusted P-value')
                ax.set_title(f'Volcano Plot: {comparison}')
                
                # Add significance lines
                ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
                ax.axvline(x=0.5, color='black', linestyle='--', alpha=0.5)
                ax.axvline(x=-0.5, color='black', linestyle='--', alpha=0.5)
                
                plt.tight_layout()
                plt.savefig(Path('results/figures') / f'volcano_{comparison}.png', 
                           dpi=300, bbox_inches='tight')
                plt.close()
    
    def _plot_de_heatmaps(self, adata, results):
        """Create heatmaps of DE genes"""
        
        if 'cell_type_markers' in results:
            # Get top marker genes for heatmap
            top_genes = []
            cell_types = []
            
            for cell_type, markers in results['cell_type_markers'].items():
                if isinstance(markers, pd.DataFrame) and len(markers) > 0:
                    top_markers = markers.head(3)['gene'].tolist()
                    top_genes.extend(top_markers)
                    cell_types.extend([cell_type] * len(top_markers))
            
            if len(top_genes) > 0:
                # Remove duplicates
                unique_genes = list(dict.fromkeys(top_genes))
                
                try:
                    # Create heatmap
                    fig, ax = plt.subplots(figsize=(12, 8))
                    
                    sc.pl.heatmap(
                        adata,
                        unique_genes[:30],  # Limit to 30 genes
                        groupby='cell_type',
                        ax=ax,
                        show=False,
                        cmap='viridis'
                    )
                    
                    plt.title('Top Marker Genes Heatmap')
                    plt.tight_layout()
                    plt.savefig(Path('results/figures') / 'marker_genes_heatmap.png', 
                               dpi=300, bbox_inches='tight')
                    plt.close()
                    
                except Exception as e:
                    self.logger.warning(f"Could not create heatmap: {e}")
                    plt.close()
