"""
Cell Annotation Module
Implements automated and manual cell type annotation methods
"""

import numpy as np
import pandas as pd
import scanpy as sc
import logging
from pathlib import Path
import matplotlib.pyplot as plt

class CellAnnotator:
    """Handles cell type annotation using various methods"""
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
    def annotate_cells(self, adata):
        """Main cell annotation function"""
        self.logger.info("Starting cell type annotation...")
        
        # Make a copy
        adata = adata.copy()
        
        # Choose annotation method
        method = self.config.get('method', 'celltypist')
        
        if method == 'celltypist':
            adata = self._celltypist_annotation(adata)
        elif method == 'manual':
            adata = self._manual_annotation(adata)
        elif method == 'scvi_label_transfer':
            adata = self._scvi_label_transfer(adata)
        elif method == 'marker_genes':
            adata = self._marker_gene_annotation(adata)
        else:
            self.logger.warning(f"Unknown annotation method: {method}, using manual")
            adata = self._manual_annotation(adata)
        
        # Validate annotations
        adata = self._validate_annotations(adata)
        
        # Plot annotation results
        self._plot_annotation_results(adata)
        
        self.logger.info("Cell annotation completed")
        
        return adata
    
    def _celltypist_annotation(self, adata):
        """Automated annotation using CellTypist"""
        try:
            import celltypist
            from celltypist import models
            
            # Ensure data is properly formatted
            adata_temp = adata.copy()
            
            # CellTypist expects raw counts normalized to 10,000 and log1p transformed
            if adata_temp.raw is not None:
                # Use raw data
                adata_temp.X = adata_temp.raw.X
                sc.pp.normalize_total(adata_temp, target_sum=1e4)
                sc.pp.log1p(adata_temp)
            else:
                # Check if data needs normalization
                max_val = adata_temp.X.max()
                if max_val > 50:  # Likely raw counts
                    sc.pp.normalize_total(adata_temp, target_sum=1e4)
                    sc.pp.log1p(adata_temp)
            
            # Make variable names unique
            adata_temp.var_names_make_unique()
            
            # Choose model based on tissue type
            tissue_type = self.config.get('tissue_type', 'general')
            
            if tissue_type in ['immune', 'blood']:
                model_name = 'Immune_All_High.pkl'
            elif tissue_type == 'lung':
                model_name = 'Adult_Lung.pkl'
            elif tissue_type == 'brain':
                model_name = 'Adult_Brain.pkl'
            else:
                model_name = 'Adult_All.pkl'  # General model
            
            # Download and load model
            models.download_models(model=[model_name])
            model = models.Model.load(model=model_name)
            
            # Predict cell types
            predictions = celltypist.annotate(adata_temp, model=model)
            
            # Transfer annotations to original data
            adata.obs['celltypist_prediction'] = predictions.predicted_labels.predicted_labels
            adata.obs['celltypist_confidence'] = predictions.predicted_labels.conf_score
            adata.obs['cell_type'] = predictions.predicted_labels.predicted_labels
            
            self.logger.info("CellTypist annotation completed")
            
        except Exception as e:
            self.logger.error(f"CellTypist annotation failed: {e}")
            # Fallback to manual annotation
            adata = self._manual_annotation(adata)
            
        return adata
    
    def _manual_annotation(self, adata):
        """Manual annotation based on cluster markers"""
        
        # Ensure we have clusters
        if 'leiden' not in adata.obs.columns:
            sc.tl.leiden(adata, resolution=0.5)
        
        # Find marker genes for each cluster
        sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
        
        # Get marker genes
        marker_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
        
        # Load cell type markers from data lake if available
        cell_type_mapping = self._load_cell_type_markers()
        
        # Initialize cell type annotations
        adata.obs['cell_type'] = 'Unknown'
        
        # Simple rule-based annotation (this would be more sophisticated in practice)
        for cluster in adata.obs['leiden'].unique():
            cluster_markers = marker_genes[cluster].head(10).tolist()
            
            # Simple matching logic (would be more complex in real implementation)
            if any(gene in ['CD3D', 'CD3E', 'CD3G'] for gene in cluster_markers):
                cell_type = 'T cell'
            elif any(gene in ['CD79A', 'CD79B', 'MS4A1'] for gene in cluster_markers):
                cell_type = 'B cell'
            elif any(gene in ['CD68', 'CD14', 'LYZ'] for gene in cluster_markers):
                cell_type = 'Macrophage'
            elif any(gene in ['CD34', 'KIT'] for gene in cluster_markers):
                cell_type = 'Stem cell'
            else:
                cell_type = f'Cluster_{cluster}'
                
            adata.obs.loc[adata.obs['leiden'] == cluster, 'cell_type'] = cell_type
        
        self.logger.info("Manual annotation completed")
        
        return adata
    
    def _scvi_label_transfer(self, adata):
        """Label transfer using scVI"""
        try:
            import scvi
            
            reference_path = self.config.get('reference_data')
            if reference_path and Path(reference_path).exists():
                # Load reference data
                adata_ref = sc.read_h5ad(reference_path)
                
                # Setup and train scANVI for label transfer
                scvi.model.SCANVI.setup_anndata(
                    adata_ref,
                    labels_key='cell_type',
                    unlabeled_category='Unknown'
                )
                
                model = scvi.model.SCANVI(adata_ref)
                model.train()
                
                # Transfer labels to query data
                adata_query = scvi.model.SCANVI.prepare_query_anndata(adata, model)
                predictions = model.predict(adata_query)
                
                adata.obs['cell_type'] = predictions
                
                self.logger.info("scVI label transfer completed")
                
            else:
                self.logger.warning("Reference data not found, using manual annotation")
                adata = self._manual_annotation(adata)
                
        except Exception as e:
            self.logger.error(f"scVI label transfer failed: {e}")
            adata = self._manual_annotation(adata)
            
        return adata
    
    def _marker_gene_annotation(self, adata):
        """Annotation using marker gene expression"""
        
        # Load marker genes from data lake
        marker_genes_dict = self._load_cell_type_markers()
        
        if not marker_genes_dict:
            self.logger.warning("No marker genes found, using manual annotation")
            return self._manual_annotation(adata)
        
        # Calculate marker gene scores
        for cell_type, markers in marker_genes_dict.items():
            # Filter markers that exist in the data
            available_markers = [gene for gene in markers if gene in adata.var_names]
            
            if len(available_markers) > 0:
                sc.tl.score_genes(adata, available_markers, score_name=f'{cell_type}_score')
        
        # Assign cell types based on highest scores
        score_columns = [col for col in adata.obs.columns if col.endswith('_score')]
        
        if score_columns:
            scores_df = adata.obs[score_columns]
            adata.obs['cell_type'] = scores_df.idxmax(axis=1).str.replace('_score', '')
        else:
            adata = self._manual_annotation(adata)
        
        return adata
    
    def _load_cell_type_markers(self):
        """Load cell type markers from data lake or config"""
        
        try:
            # Try to load from data lake
            import boto3
            import pandas as pd
            
            s3 = boto3.client('s3')
            
            # Download marker file
            s3.download_file(
                'biomni-datalake', 
                'marker_celltype.parquet',
                '/tmp/marker_celltype.parquet'
            )
            
            markers_df = pd.read_parquet('/tmp/marker_celltype.parquet')
            
            # Convert to dictionary format
            marker_dict = {}
            for cell_type in markers_df['cell_type'].unique():
                markers = markers_df[markers_df['cell_type'] == cell_type]['gene'].tolist()
                marker_dict[cell_type] = markers
            
            return marker_dict
            
        except Exception as e:
            self.logger.warning(f"Could not load markers from data lake: {e}")
            
            # Use default marker genes
            default_markers = {
                'T cell': ['CD3D', 'CD3E', 'CD3G', 'CD2', 'CD8A', 'CD4'],
                'B cell': ['CD79A', 'CD79B', 'MS4A1', 'CD19', 'IGHM'],
                'NK cell': ['GNLY', 'NKG7', 'CD160', 'KLRD1'],
                'Macrophage': ['CD68', 'CD14', 'LYZ', 'CD163', 'MRC1'],
                'Dendritic cell': ['CD1C', 'FCER1A', 'CLEC10A', 'CD1E'],
                'Monocyte': ['CD14', 'LYZ', 'S100A9', 'FCN1'],
                'Neutrophil': ['FCGR3B', 'CSF3R', 'FPR1'],
                'Endothelial cell': ['PECAM1', 'VWF', 'ENG', 'CDH5'],
                'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM'],
                'Epithelial cell': ['EPCAM', 'KRT8', 'KRT18', 'KRT19']
            }
            
            return default_markers
    
    def _validate_annotations(self, adata):
        """Validate cell type annotations"""
        
        if 'cell_type' not in adata.obs.columns:
            self.logger.error("No cell type annotations found")
            return adata
        
        # Check for unknown/unassigned cells
        unknown_cells = (adata.obs['cell_type'] == 'Unknown').sum()
        total_cells = adata.shape[0]
        unknown_percentage = unknown_cells / total_cells * 100
        
        self.logger.info(f"Annotation coverage: {100-unknown_percentage:.1f}% ({total_cells-unknown_cells}/{total_cells} cells)")
        
        if unknown_percentage > 50:
            self.logger.warning(f"High percentage of unknown cells: {unknown_percentage:.1f}%")
        
        return adata
    
    def _plot_annotation_results(self, adata):
        """Plot annotation results"""
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # UMAP colored by cell type
        sc.pl.umap(adata, color='cell_type', ax=axes[0, 0], show=False)
        axes[0, 0].set_title('Cell Types')
        
        # UMAP colored by clusters
        if 'leiden' in adata.obs.columns:
            sc.pl.umap(adata, color='leiden', ax=axes[0, 1], show=False, legend_loc='right margin')
            axes[0, 1].set_title('Leiden Clusters')
        
        # Cell type counts
        cell_type_counts = adata.obs['cell_type'].value_counts()
        axes[1, 0].barh(range(len(cell_type_counts)), cell_type_counts.values)
        axes[1, 0].set_yticks(range(len(cell_type_counts)))
        axes[1, 0].set_yticklabels(cell_type_counts.index)
        axes[1, 0].set_xlabel('Number of cells')
        axes[1, 0].set_title('Cell type distribution')
        
        # Cluster-cell type confusion matrix
        if 'leiden' in adata.obs.columns:
            confusion_matrix = pd.crosstab(adata.obs['leiden'], adata.obs['cell_type'])
            im = axes[1, 1].imshow(confusion_matrix.values, aspect='auto', cmap='Blues')
            axes[1, 1].set_xticks(range(len(confusion_matrix.columns)))
            axes[1, 1].set_xticklabels(confusion_matrix.columns, rotation=45, ha='right')
            axes[1, 1].set_yticks(range(len(confusion_matrix.index)))
            axes[1, 1].set_yticklabels(confusion_matrix.index)
            axes[1, 1].set_xlabel('Cell Type')
            axes[1, 1].set_ylabel('Leiden Cluster')
            axes[1, 1].set_title('Cluster-Cell Type Matrix')
            plt.colorbar(im, ax=axes[1, 1])
        
        plt.tight_layout()
        
        # Save plot
        output_dir = Path('results/figures')
        output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_dir / 'cell_annotation_results.png', dpi=300, bbox_inches='tight')
        plt.close()
