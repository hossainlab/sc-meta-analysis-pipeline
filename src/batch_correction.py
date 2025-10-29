"""
Batch Correction Module
Implements various batch correction methods for single-cell data integration
"""

import numpy as np
import pandas as pd
import scanpy as sc
import logging
from pathlib import Path
import matplotlib.pyplot as plt

class BatchCorrector:
    """Handles batch correction and data integration"""
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
    def integrate_datasets(self, datasets):
        """Main integration function"""
        self.logger.info(f"Integrating {len(datasets)} datasets...")
        
        # Concatenate datasets
        adata_concat = self._concatenate_datasets(datasets)
        
        # Perform batch correction
        adata_integrated = self._correct_batch_effects(adata_concat)
        
        # Generate integration plots
        self._plot_integration_results(adata_integrated)
        
        self.logger.info("Integration completed")
        
        return adata_integrated
    
    def _concatenate_datasets(self, datasets):
        """Concatenate multiple datasets"""
        
        adata_list = list(datasets.values())
        dataset_names = list(datasets.keys())
        
        # Ensure all datasets have the same genes
        common_genes = set(adata_list[0].var_names)
        for adata in adata_list[1:]:
            common_genes = common_genes.intersection(set(adata.var_names))
        
        self.logger.info(f"Found {len(common_genes)} common genes across datasets")
        
        # Subset to common genes
        for i in range(len(adata_list)):
            adata_list[i] = adata_list[i][:, list(common_genes)].copy()
        
        # Concatenate
        adata_concat = adata_list[0].concatenate(
            adata_list[1:],
            batch_key='dataset',
            batch_categories=dataset_names
        )
        
        return adata_concat
    
    def _correct_batch_effects(self, adata):
        """Apply batch correction methods"""
        
        method = self.config.get('method', 'harmony')
        
        if method == 'harmony':
            adata_corrected = self._harmony_integration(adata)
        elif method == 'scanorama':
            adata_corrected = self._scanorama_integration(adata)
        elif method == 'scvi':
            adata_corrected = self._scvi_integration(adata)
        elif method == 'combat':
            adata_corrected = self._combat_integration(adata)
        else:
            self.logger.warning(f"Unknown batch correction method: {method}")
            adata_corrected = adata
            
        return adata_corrected
    
    def _harmony_integration(self, adata):
        """Harmony batch correction"""
        try:
            # Compute PCA first
            sc.tl.pca(adata, svd_solver='arpack')
            
            # Try to use harmony-pytorch
            try:
                import harmony
                harmony_out = harmony.run_harmony(
                    adata.obsm['X_pca'],
                    adata.obs,
                    ['dataset'],
                    max_iter_harmony=20
                )
                adata.obsm['X_pca_harmony'] = harmony_out.Z_corr.T
                
            except ImportError:
                self.logger.warning("harmony-pytorch not available, using scanpy's external harmony")
                sc.external.pp.harmony_integrate(
                    adata, 
                    key='dataset',
                    basis='X_pca',
                    adjusted_basis='X_pca_harmony'
                )
            
            # Compute neighbors and UMAP on corrected data
            sc.pp.neighbors(adata, use_rep='X_pca_harmony')
            sc.tl.umap(adata)
            sc.tl.leiden(adata, resolution=0.5)
            
            self.logger.info("Harmony integration completed")
            
        except Exception as e:
            self.logger.error(f"Harmony integration failed: {e}")
            # Fallback to standard processing
            sc.tl.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            sc.tl.leiden(adata, resolution=0.5)
            
        return adata
    
    def _scanorama_integration(self, adata):
        """Scanorama batch correction"""
        try:
            import scanorama
            
            # Split by batches
            batches = []
            batch_names = []
            
            for batch in adata.obs['dataset'].unique():
                batch_data = adata[adata.obs['dataset'] == batch]
                batches.append(batch_data.X.toarray() if hasattr(batch_data.X, 'toarray') else batch_data.X)
                batch_names.append(str(batch))
            
            # Run scanorama
            integrated, corrected_genes = scanorama.integrate(
                batches,
                [adata.var_names.tolist()] * len(batches)
            )
            
            # Reconstruct integrated data
            X_integrated = np.vstack(integrated)
            adata.X = X_integrated
            
            # Compute PCA, neighbors, UMAP
            sc.tl.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            sc.tl.leiden(adata, resolution=0.5)
            
            self.logger.info("Scanorama integration completed")
            
        except Exception as e:
            self.logger.error(f"Scanorama integration failed: {e}")
            # Fallback
            sc.tl.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            sc.tl.leiden(adata, resolution=0.5)
            
        return adata
    
    def _scvi_integration(self, adata):
        """scVI-tools integration"""
        try:
            import scvi
            
            # Setup scVI model
            scvi.model.SCVI.setup_anndata(
                adata, 
                layer=None,
                batch_key='dataset'
            )
            
            model = scvi.model.SCVI(adata)
            model.train()
            
            # Get latent representation
            adata.obsm["X_scVI"] = model.get_latent_representation()
            
            # Compute neighbors and UMAP
            sc.pp.neighbors(adata, use_rep="X_scVI")
            sc.tl.umap(adata)
            sc.tl.leiden(adata, resolution=0.5)
            
            self.logger.info("scVI integration completed")
            
        except Exception as e:
            self.logger.error(f"scVI integration failed: {e}")
            # Fallback
            sc.tl.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            sc.tl.leiden(adata, resolution=0.5)
            
        return adata
    
    def _combat_integration(self, adata):
        """ComBat batch correction"""
        try:
            # ComBat correction
            sc.pp.combat(adata, key='dataset')
            
            # Standard processing after correction
            sc.tl.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            sc.tl.leiden(adata, resolution=0.5)
            
            self.logger.info("ComBat integration completed")
            
        except Exception as e:
            self.logger.error(f"ComBat integration failed: {e}")
            # Fallback
            sc.tl.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            sc.tl.leiden(adata, resolution=0.5)
            
        return adata
    
    def _plot_integration_results(self, adata):
        """Plot integration results"""
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # UMAP colored by dataset
        sc.pl.umap(adata, color='dataset', ax=axes[0, 0], show=False)
        axes[0, 0].set_title('Datasets')
        
        # UMAP colored by clusters
        if 'leiden' in adata.obs.columns:
            sc.pl.umap(adata, color='leiden', ax=axes[0, 1], show=False)
            axes[0, 1].set_title('Clusters')
        
        # Dataset composition per cluster
        if 'leiden' in adata.obs.columns:
            cluster_dataset_counts = pd.crosstab(adata.obs['leiden'], adata.obs['dataset'])
            cluster_dataset_proportions = cluster_dataset_counts.div(cluster_dataset_counts.sum(axis=1), axis=0)
            
            im = axes[1, 0].imshow(cluster_dataset_proportions.values, aspect='auto', cmap='Blues')
            axes[1, 0].set_xticks(range(len(cluster_dataset_proportions.columns)))
            axes[1, 0].set_xticklabels(cluster_dataset_proportions.columns, rotation=45)
            axes[1, 0].set_yticks(range(len(cluster_dataset_proportions.index)))
            axes[1, 0].set_yticklabels(cluster_dataset_proportions.index)
            axes[1, 0].set_xlabel('Dataset')
            axes[1, 0].set_ylabel('Cluster')
            axes[1, 0].set_title('Dataset composition per cluster')
            plt.colorbar(im, ax=axes[1, 0])
        
        # Cell counts per dataset
        dataset_counts = adata.obs['dataset'].value_counts()
        axes[1, 1].bar(dataset_counts.index, dataset_counts.values)
        axes[1, 1].set_xlabel('Dataset')
        axes[1, 1].set_ylabel('Number of cells')
        axes[1, 1].set_title('Cell counts per dataset')
        axes[1, 1].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        
        # Save plot
        output_dir = Path('results/figures')
        output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_dir / 'integration_results.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def calculate_integration_metrics(self, adata):
        """Calculate integration quality metrics"""
        
        metrics = {}
        
        # Silhouette score for batch mixing
        try:
            from sklearn.metrics import silhouette_score
            
            if 'X_pca' in adata.obsm:
                silhouette_batch = silhouette_score(adata.obsm['X_pca'], adata.obs['dataset'])
                metrics['silhouette_batch'] = silhouette_batch
                
            if 'leiden' in adata.obs:
                silhouette_cluster = silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'])
                metrics['silhouette_cluster'] = silhouette_cluster
                
        except Exception as e:
            self.logger.warning(f"Could not calculate silhouette scores: {e}")
        
        # Save metrics
        output_dir = Path('results/tables')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        pd.DataFrame([metrics]).to_csv(output_dir / 'integration_metrics.csv', index=False)
        
        return metrics
