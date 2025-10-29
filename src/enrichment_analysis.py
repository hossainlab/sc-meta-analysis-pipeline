"""
Enrichment Analysis Module
Implements pathway and gene set enrichment analysis
"""

import numpy as np
import pandas as pd
import logging
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

class EnrichmentAnalyzer:
    """Handles pathway enrichment and gene set analysis"""
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
    def run_gsea(self, de_results, adata):
        """Main GSEA analysis function"""
        self.logger.info("Starting enrichment analysis...")
        
        enrichment_results = {}
        
        # Load gene sets
        gene_sets = self._load_gene_sets()
        
        # Run GSEA for different comparisons
        for analysis_type, results in de_results.items():
            if isinstance(results, dict):
                for comparison, de_df in results.items():
                    if isinstance(de_df, pd.DataFrame):
                        gsea_result = self._run_gsea_analysis(de_df, gene_sets, f"{analysis_type}_{comparison}")
                        if gsea_result is not None:
                            enrichment_results[f"{analysis_type}_{comparison}"] = gsea_result
        
        # Run AUCell for cell-level gene set activity
        aucell_results = self._run_aucell_analysis(adata, gene_sets)
        if aucell_results is not None:
            enrichment_results['aucell'] = aucell_results
        
        # Save results
        self._save_enrichment_results(enrichment_results)
        
        # Create visualizations
        self._plot_enrichment_results(enrichment_results, adata)
        
        self.logger.info("Enrichment analysis completed")
        
        return enrichment_results
    
    def _load_gene_sets(self):
        """Load gene sets from various sources"""
        
        gene_sets = {}
        
        # Try to load from decoupler
        try:
            import decoupler as dc
            
            # MSigDB gene sets
            msigdb = dc.get_resource('MSigDB')
            gene_sets['msigdb'] = msigdb
            
        except Exception as e:
            self.logger.warning(f"Could not load MSigDB from decoupler: {e}")
        
        # Try to load from data lake
        try:
            import boto3
            import pandas as pd
            
            s3 = boto3.client('s3')
            
            # Load oncogenic signatures
            s3.download_file(
                'biomni-datalake',
                'msigdb_human_c6_oncogenic_signature_geneset.parquet',
                '/tmp/oncogenic.parquet'
            )
            oncogenic_df = pd.read_parquet('/tmp/oncogenic.parquet')
            gene_sets['oncogenic'] = oncogenic_df
            
            # Load immunologic signatures  
            s3.download_file(
                'biomni-datalake',
                'msigdb_human_c7_immunologic_signature_geneset.parquet',
                '/tmp/immunologic.parquet'
            )
            immuno_df = pd.read_parquet('/tmp/immunologic.parquet')
            gene_sets['immunologic'] = immuno_df
            
            # Load cell type signatures
            s3.download_file(
                'biomni-datalake',
                'msigdb_human_c8_celltype_signature_geneset.parquet',
                '/tmp/celltype.parquet'
            )
            celltype_df = pd.read_parquet('/tmp/celltype.parquet')
            gene_sets['celltype'] = celltype_df
            
        except Exception as e:
            self.logger.warning(f"Could not load gene sets from data lake: {e}")
        
        # Default gene sets if nothing else works
        if not gene_sets:
            gene_sets = self._get_default_gene_sets()
        
        return gene_sets
    
    def _get_default_gene_sets(self):
        """Get default gene sets"""
        
        default_sets = {
            'default': pd.DataFrame({
                'source': ['Immune', 'Immune', 'Immune', 'Cell_Cycle', 'Cell_Cycle'],
                'target': ['CD3D', 'CD3E', 'CD3G', 'MKI67', 'PCNA']
            })
        }
        
        return default_sets
    
    def _run_gsea_analysis(self, de_df, gene_sets, comparison_name):
        """Run GSEA analysis using decoupler"""
        
        try:
            import decoupler as dc
            
            # Prepare gene statistics (t-statistics or log fold changes)
            if 'logfoldchanges' in de_df.columns:
                gene_stats = de_df.set_index('gene')['logfoldchanges'].to_dict()
            elif 'scores' in de_df.columns:
                gene_stats = de_df.set_index('gene')['scores'].to_dict()
            else:
                self.logger.warning(f"No suitable statistics column found for {comparison_name}")
                return None
            
            # Convert to DataFrame format expected by decoupler
            stats_df = pd.DataFrame(list(gene_stats.items()), columns=['gene', 'stat'])
            stats_df = stats_df.set_index('gene').T
            
            gsea_results = {}
            
            for gene_set_name, gene_set in gene_sets.items():
                if isinstance(gene_set, pd.DataFrame):
                    
                    # Ensure proper column names
                    if 'source' in gene_set.columns and 'target' in gene_set.columns:
                        net = gene_set.copy()
                    else:
                        # Try to infer column structure
                        columns = gene_set.columns.tolist()
                        if len(columns) >= 2:
                            net = gene_set.copy()
                            net.columns = ['source', 'target'] + columns[2:]
                        else:
                            continue
                    
                    # Filter gene sets by size
                    set_sizes = net.groupby('source').size()
                    valid_sets = set_sizes[(set_sizes >= 15) & (set_sizes <= 500)].index
                    net = net[net['source'].isin(valid_sets)]
                    
                    if len(net) > 0:
                        try:
                            # Run GSEA
                            gsea_result = dc.run_gsea(
                                mat=stats_df,
                                net=net,
                                source='source',
                                target='target',
                                min_n=15,
                                verbose=False
                            )
                            
                            gsea_results[gene_set_name] = gsea_result
                            
                        except Exception as e:
                            self.logger.warning(f"GSEA failed for {gene_set_name}: {e}")
            
            return gsea_results if gsea_results else None
            
        except Exception as e:
            self.logger.error(f"GSEA analysis failed for {comparison_name}: {e}")
            return None
    
    def _run_aucell_analysis(self, adata, gene_sets):
        """Run AUCell analysis for single-cell gene set activity"""
        
        try:
            import decoupler as dc
            
            aucell_results = {}
            
            for gene_set_name, gene_set in gene_sets.items():
                if isinstance(gene_set, pd.DataFrame):
                    
                    # Ensure proper column names
                    if 'source' in gene_set.columns and 'target' in gene_set.columns:
                        net = gene_set.copy()
                    else:
                        columns = gene_set.columns.tolist()
                        if len(columns) >= 2:
                            net = gene_set.copy()
                            net.columns = ['source', 'target'] + columns[2:]
                        else:
                            continue
                    
                    # Filter gene sets
                    set_sizes = net.groupby('source').size()
                    valid_sets = set_sizes[(set_sizes >= 15) & (set_sizes <= 500)].index
                    net = net[net['source'].isin(valid_sets)]
                    
                    if len(net) > 0:
                        try:
                            # Run AUCell
                            auc_result = dc.run_aucell(
                                mat=adata,
                                net=net,
                                source='source',
                                target='target',
                                verbose=False
                            )
                            
                            aucell_results[gene_set_name] = auc_result
                            
                        except Exception as e:
                            self.logger.warning(f"AUCell failed for {gene_set_name}: {e}")
            
            return aucell_results if aucell_results else None
            
        except Exception as e:
            self.logger.error(f"AUCell analysis failed: {e}")
            return None
    
    def _save_enrichment_results(self, results):
        """Save enrichment analysis results"""
        
        output_dir = Path('results/tables')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        for analysis_name, analysis_results in results.items():
            if isinstance(analysis_results, dict):
                for gene_set_name, result in analysis_results.items():
                    if hasattr(result, 'to_csv'):
                        filename = f"enrichment_{analysis_name}_{gene_set_name}.csv"
                        result.to_csv(output_dir / filename)
                    elif isinstance(result, tuple) and len(result) >= 3:
                        # Decoupler results format (estimates, pvals, norms)
                        estimates, pvals, norms = result[:3]
                        if hasattr(estimates, 'to_csv'):
                            estimates.to_csv(output_dir / f"enrichment_{analysis_name}_{gene_set_name}_estimates.csv")
                        if hasattr(pvals, 'to_csv'):
                            pvals.to_csv(output_dir / f"enrichment_{analysis_name}_{gene_set_name}_pvals.csv")
        
        self.logger.info("Enrichment results saved")
    
    def _plot_enrichment_results(self, results, adata):
        """Create enrichment analysis visualizations"""
        
        output_dir = Path('results/figures')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Plot GSEA results
        self._plot_gsea_results(results)
        
        # Plot AUCell results
        if 'aucell' in results:
            self._plot_aucell_results(results['aucell'], adata)
    
    def _plot_gsea_results(self, results):
        """Plot GSEA results"""
        
        for analysis_name, analysis_results in results.items():
            if 'aucell' in analysis_name:
                continue
                
            if isinstance(analysis_results, dict):
                for gene_set_name, result in list(analysis_results.items())[:2]:  # Limit plots
                    
                    try:
                        if isinstance(result, tuple) and len(result) >= 2:
                            estimates, pvals = result[:2]
                            
                            if hasattr(estimates, 'iloc') and hasattr(pvals, 'iloc'):
                                # Create a combined dataframe for plotting
                                plot_df = pd.DataFrame({
                                    'pathway': estimates.index,
                                    'nes': estimates.iloc[:, 0] if estimates.shape[1] > 0 else 0,
                                    'pval': pvals.iloc[:, 0] if pvals.shape[1] > 0 else 1
                                })
                                
                                # Filter significant pathways
                                sig_pathways = plot_df[plot_df['pval'] < 0.05].head(20)
                                
                                if len(sig_pathways) > 0:
                                    fig, ax = plt.subplots(figsize=(10, 8))
                                    
                                    # Create bar plot
                                    bars = ax.barh(range(len(sig_pathways)), sig_pathways['nes'])
                                    
                                    # Color bars by NES direction
                                    for i, (bar, nes) in enumerate(zip(bars, sig_pathways['nes'])):
                                        bar.set_color('red' if nes > 0 else 'blue')
                                    
                                    ax.set_yticks(range(len(sig_pathways)))
                                    ax.set_yticklabels([p[:50] for p in sig_pathways['pathway']], fontsize=8)
                                    ax.set_xlabel('Normalized Enrichment Score')
                                    ax.set_title(f'GSEA Results: {analysis_name} - {gene_set_name}')
                                    ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
                                    
                                    plt.tight_layout()
                                    plt.savefig(
                                        output_dir / f'gsea_{analysis_name}_{gene_set_name}.png',
                                        dpi=300, bbox_inches='tight'
                                    )
                                    plt.close()
                    
                    except Exception as e:
                        self.logger.warning(f"Could not plot GSEA results for {analysis_name}-{gene_set_name}: {e}")
    
    def _plot_aucell_results(self, aucell_results, adata):
        """Plot AUCell results"""
        
        try:
            for gene_set_name, result in list(aucell_results.items())[:2]:  # Limit plots
                if isinstance(result, tuple) and len(result) >= 1:
                    auc_scores = result[0]
                    
                    if hasattr(auc_scores, 'T'):
                        # Add AUC scores to adata for plotting
                        for pathway in auc_scores.columns[:5]:  # Plot top 5 pathways
                            adata.obs[f'AUC_{pathway}'] = auc_scores[pathway].values
                    
                        # Create UMAP plots
                        import scanpy as sc
                        
                        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
                        axes = axes.flatten()
                        
                        for i, pathway in enumerate(auc_scores.columns[:5]):
                            if i < len(axes):
                                try:
                                    sc.pl.umap(
                                        adata,
                                        color=f'AUC_{pathway}',
                                        ax=axes[i],
                                        show=False,
                                        cmap='viridis'
                                    )
                                    axes[i].set_title(pathway[:30])
                                except:
                                    axes[i].text(0.5, 0.5, 'Plot failed', ha='center', va='center')
                        
                        # Remove extra subplots
                        for j in range(i+1, len(axes)):
                            fig.delaxes(axes[j])
                        
                        plt.suptitle(f'AUCell Results: {gene_set_name}')
                        plt.tight_layout()
                        plt.savefig(
                            Path('results/figures') / f'aucell_{gene_set_name}.png',
                            dpi=300, bbox_inches='tight'
                        )
                        plt.close()
                        
        except Exception as e:
            self.logger.warning(f"Could not plot AUCell results: {e}")
