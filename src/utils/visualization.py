"""
Visualization utilities for the pipeline
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path
from typing import Optional, List, Dict, Union
import logging

def setup_plotting_style(style: str = 'default', dpi: int = 300):
    """Setup consistent plotting style"""
    
    if style == 'publication':
        plt.style.use('seaborn-v0_8-paper')
        params = {
            'figure.figsize': (8, 6),
            'font.size': 12,
            'axes.titlesize': 14,
            'axes.labelsize': 12,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 10,
            'figure.dpi': dpi,
            'savefig.dpi': dpi,
            'savefig.bbox': 'tight'
        }
        plt.rcParams.update(params)
    else:
        sc.settings.set_figure_params(dpi=dpi, facecolor='white')

def create_qc_dashboard(adata: sc.AnnData, 
                       output_dir: Union[str, Path],
                       sample_name: str = 'sample') -> Path:
    """Create comprehensive QC dashboard"""
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    fig.suptitle(f'Quality Control Dashboard - {sample_name}', fontsize=16)
    
    # 1. Number of genes per cell
    adata.obs['n_genes'].hist(bins=50, ax=axes[0, 0])
    axes[0, 0].set_xlabel('Number of genes')
    axes[0, 0].set_ylabel('Number of cells')
    axes[0, 0].set_title('Genes per cell distribution')
    
    # 2. Total counts per cell
    adata.obs['total_counts'].hist(bins=50, ax=axes[0, 1])
    axes[0, 1].set_xlabel('Total counts')
    axes[0, 1].set_ylabel('Number of cells')
    axes[0, 1].set_title('UMI counts per cell distribution')
    axes[0, 1].set_xscale('log')
    
    # 3. Mitochondrial gene percentage
    if 'pct_counts_mt' in adata.obs.columns:
        adata.obs['pct_counts_mt'].hist(bins=50, ax=axes[0, 2])
        axes[0, 2].set_xlabel('Mitochondrial gene %')
        axes[0, 2].set_ylabel('Number of cells')
        axes[0, 2].set_title('Mitochondrial gene percentage')
    
    # 4. Genes vs UMI scatter
    axes[1, 0].scatter(adata.obs['total_counts'], adata.obs['n_genes'], 
                      alpha=0.5, s=1)
    axes[1, 0].set_xlabel('Total counts')
    axes[1, 0].set_ylabel('Number of genes')
    axes[1, 0].set_title('Genes vs UMI counts')
    axes[1, 0].set_xscale('log')
    
    # 5. UMI vs Mitochondrial percentage
    if 'pct_counts_mt' in adata.obs.columns:
        axes[1, 1].scatter(adata.obs['total_counts'], adata.obs['pct_counts_mt'], 
                          alpha=0.5, s=1)
        axes[1, 1].set_xlabel('Total counts')
        axes[1, 1].set_ylabel('Mitochondrial gene %')
        axes[1, 1].set_title('UMI counts vs Mitochondrial %')
        axes[1, 1].set_xscale('log')
    
    # 6. Gene detection statistics
    gene_detection = (adata.X > 0).sum(axis=0)
    if hasattr(gene_detection, 'A1'):
        gene_detection = gene_detection.A1
    
    axes[1, 2].hist(gene_detection, bins=50)
    axes[1, 2].set_xlabel('Number of cells')
    axes[1, 2].set_ylabel('Number of genes')
    axes[1, 2].set_title('Gene detection across cells')
    axes[1, 2].set_xscale('log')
    
    # 7. Top expressed genes
    top_genes = pd.Series(adata.X.sum(axis=0)).sort_values(ascending=False)
    if hasattr(top_genes, 'A1'):
        top_genes = pd.Series(top_genes.A1, index=adata.var_names).sort_values(ascending=False)
    
    top_genes.head(10).plot(kind='barh', ax=axes[2, 0])
    axes[2, 0].set_xlabel('Total expression')
    axes[2, 0].set_title('Top 10 expressed genes')
    
    # 8. Cell type distribution (if available)
    if 'cell_type' in adata.obs.columns:
        cell_type_counts = adata.obs['cell_type'].value_counts()
        cell_type_counts.plot(kind='pie', ax=axes[2, 1], autopct='%1.1f%%')
        axes[2, 1].set_title('Cell type distribution')
        axes[2, 1].set_ylabel('')
    else:
        axes[2, 1].text(0.5, 0.5, 'Cell types not annotated', 
                       ha='center', va='center', transform=axes[2, 1].transAxes)
        axes[2, 1].set_title('Cell type distribution')
    
    # 9. Batch distribution (if available)
    if 'batch' in adata.obs.columns or 'dataset' in adata.obs.columns:
        batch_col = 'batch' if 'batch' in adata.obs.columns else 'dataset'
        batch_counts = adata.obs[batch_col].value_counts()
        batch_counts.plot(kind='bar', ax=axes[2, 2])
        axes[2, 2].set_xlabel('Batch/Dataset')
        axes[2, 2].set_ylabel('Number of cells')
        axes[2, 2].set_title('Batch distribution')
        axes[2, 2].tick_params(axis='x', rotation=45)
    else:
        axes[2, 2].text(0.5, 0.5, 'No batch information', 
                       ha='center', va='center', transform=axes[2, 2].transAxes)
        axes[2, 2].set_title('Batch distribution')
    
    plt.tight_layout()
    
    # Save dashboard
    dashboard_path = output_path / f'{sample_name}_qc_dashboard.png'
    plt.savefig(dashboard_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return dashboard_path

def create_integration_plot(adata: sc.AnnData,
                           output_dir: Union[str, Path],
                           color_by: List[str] = None) -> List[Path]:
    """Create integration visualization plots"""
    
    if color_by is None:
        color_by = ['dataset', 'cell_type', 'leiden']
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    saved_plots = []
    
    # Filter color_by to only include available columns
    available_colors = [col for col in color_by if col in adata.obs.columns]
    
    if len(available_colors) == 0:
        logging.warning("No valid columns found for coloring integration plots")
        return saved_plots
    
    # Create subplot grid
    n_plots = len(available_colors)
    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
    
    if n_plots == 1:
        axes = [axes]
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    
    for i, color in enumerate(available_colors):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col] if n_rows > 1 else axes[col]
        
        sc.pl.umap(adata, color=color, ax=ax, show=False, frameon=False)
        ax.set_title(f'Colored by {color}')
    
    # Hide extra subplots
    for i in range(n_plots, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row, col].set_visible(False)
    
    plt.tight_layout()
    
    plot_path = output_path / 'integration_overview.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    saved_plots.append(plot_path)
    
    return saved_plots

def create_marker_gene_plots(adata: sc.AnnData,
                           marker_genes: Dict[str, List[str]],
                           output_dir: Union[str, Path]) -> List[Path]:
    """Create marker gene visualization plots"""
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    saved_plots = []
    
    # Flatten marker genes and check availability
    all_markers = []
    for cell_type, genes in marker_genes.items():
        all_markers.extend(genes)
    
    # Find available markers
    available_markers = [gene for gene in all_markers if gene in adata.var_names]
    
    if len(available_markers) == 0:
        logging.warning("No marker genes found in dataset")
        return saved_plots
    
    # Create dotplot
    if 'cell_type' in adata.obs.columns:
        try:
            fig, ax = plt.subplots(figsize=(12, 8))
            sc.pl.dotplot(adata, available_markers[:30], groupby='cell_type', 
                         ax=ax, show=False)
            plt.title('Marker Gene Expression')
            plt.tight_layout()
            
            dotplot_path = output_path / 'marker_genes_dotplot.png'
            plt.savefig(dotplot_path, dpi=300, bbox_inches='tight')
            plt.close()
            saved_plots.append(dotplot_path)
            
        except Exception as e:
            logging.warning(f"Could not create dotplot: {e}")
    
    # Create UMAP plots for top markers
    top_markers = available_markers[:6]  # Show top 6 markers
    
    if len(top_markers) > 0:
        n_cols = 3
        n_rows = (len(top_markers) + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
        
        if len(top_markers) == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.reshape(1, -1)
        
        for i, gene in enumerate(top_markers):
            row = i // n_cols
            col = i % n_cols
            ax = axes[row, col] if n_rows > 1 else axes[col]
            
            sc.pl.umap(adata, color=gene, ax=ax, show=False, frameon=False)
            ax.set_title(f'{gene} expression')
        
        # Hide extra subplots
        for i in range(len(top_markers), n_rows * n_cols):
            row = i // n_cols
            col = i % n_cols
            axes[row, col].set_visible(False)
        
        plt.tight_layout()
        
        umap_path = output_path / 'marker_genes_umap.png'
        plt.savefig(umap_path, dpi=300, bbox_inches='tight')
        plt.close()
        saved_plots.append(umap_path)
    
    return saved_plots

def save_publication_plots(adata: sc.AnnData, 
                         output_dir: Union[str, Path],
                         plot_types: List[str] = None) -> Dict[str, Path]:
    """Create publication-ready plots"""
    
    if plot_types is None:
        plot_types = ['umap_celltype', 'umap_dataset', 'qc_metrics']
    
    output_path = Path(output_dir) / 'publication'
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Setup publication style
    setup_plotting_style('publication')
    
    saved_plots = {}
    
    for plot_type in plot_types:
        try:
            if plot_type == 'umap_celltype' and 'cell_type' in adata.obs.columns:
                fig, ax = plt.subplots(figsize=(8, 6))
                sc.pl.umap(adata, color='cell_type', ax=ax, show=False, 
                          frameon=False, legend_loc='right margin')
                plt.title('Cell Type Annotation', fontsize=14, pad=20)
                
                plot_path = output_path / 'umap_cell_types.pdf'
                plt.savefig(plot_path, format='pdf', bbox_inches='tight')
                plt.savefig(plot_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
                plt.close()
                
                saved_plots['umap_celltype'] = plot_path
                
            elif plot_type == 'umap_dataset' and 'dataset' in adata.obs.columns:
                fig, ax = plt.subplots(figsize=(8, 6))
                sc.pl.umap(adata, color='dataset', ax=ax, show=False, 
                          frameon=False, legend_loc='right margin')
                plt.title('Dataset Integration', fontsize=14, pad=20)
                
                plot_path = output_path / 'umap_datasets.pdf'
                plt.savefig(plot_path, format='pdf', bbox_inches='tight')
                plt.savefig(plot_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
                plt.close()
                
                saved_plots['umap_dataset'] = plot_path
                
            elif plot_type == 'qc_metrics':
                fig, axes = plt.subplots(1, 3, figsize=(15, 4))
                
                # QC violin plots
                qc_metrics = ['n_genes_by_counts', 'total_counts']
                if 'pct_counts_mt' in adata.obs.columns:
                    qc_metrics.append('pct_counts_mt')
                
                for i, metric in enumerate(qc_metrics[:3]):
                    if metric in adata.obs.columns:
                        if 'cell_type' in adata.obs.columns:
                            sns.violinplot(data=adata.obs, x='cell_type', y=metric, 
                                         ax=axes[i], inner='quartile')
                            axes[i].tick_params(axis='x', rotation=45)
                        else:
                            adata.obs[metric].hist(bins=50, ax=axes[i])
                        
                        axes[i].set_title(metric.replace('_', ' ').title())
                
                plt.tight_layout()
                
                plot_path = output_path / 'qc_metrics.pdf'
                plt.savefig(plot_path, format='pdf', bbox_inches='tight')
                plt.savefig(plot_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
                plt.close()
                
                saved_plots['qc_metrics'] = plot_path
                
        except Exception as e:
            logging.warning(f"Could not create {plot_type} plot: {e}")
    
    return saved_plots
