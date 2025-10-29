"""
Meta-Analysis Module
Combines results across multiple studies for meta-analysis
"""

import numpy as np
import pandas as pd
import logging
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import StandardScaler

class MetaAnalyzer:
    """Handles meta-analysis across multiple studies"""
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
    def combine_studies(self, datasets, de_results):
        """Main meta-analysis function"""
        self.logger.info("Starting meta-analysis...")
        
        meta_results = {}
        
        # Meta-analysis of differential expression results
        if de_results:
            meta_results['de_meta'] = self._meta_analyze_de(de_results)
        
        # Meta-analysis of cell type compositions
        meta_results['composition_meta'] = self._meta_analyze_composition(datasets)
        
        # Meta-analysis of gene expression patterns
        meta_results['expression_meta'] = self._meta_analyze_expression(datasets)
        
        # Calculate study heterogeneity
        meta_results['heterogeneity'] = self._calculate_heterogeneity(datasets)
        
        # Save results
        self._save_meta_results(meta_results)
        
        # Create visualizations
        self._plot_meta_results(meta_results, datasets)
        
        self.logger.info("Meta-analysis completed")
        
        return meta_results
    
    def _meta_analyze_de(self, de_results):
        """Meta-analysis of differential expression results"""
        
        meta_de_results = {}
        
        # Collect all DE comparisons
        all_comparisons = {}
        for study, study_results in de_results.items():
            if isinstance(study_results, dict):
                for comparison, de_df in study_results.items():
                    if isinstance(de_df, pd.DataFrame):
                        comparison_type = comparison.split('_')[-1] if '_' in comparison else comparison
                        
                        if comparison_type not in all_comparisons:
                            all_comparisons[comparison_type] = []
                        
                        # Add study identifier
                        de_df_copy = de_df.copy()
                        de_df_copy['study'] = study
                        all_comparisons[comparison_type].append(de_df_copy)
        
        # Perform meta-analysis for each comparison type
        for comp_type, comp_results in all_comparisons.items():
            if len(comp_results) > 1:  # Need at least 2 studies
                meta_result = self._combine_de_studies(comp_results, comp_type)
                if meta_result is not None:
                    meta_de_results[comp_type] = meta_result
        
        return meta_de_results
    
    def _combine_de_studies(self, study_results, comparison_type):
        """Combine DE results across studies using meta-analysis"""
        
        try:
            # Get common genes across all studies
            all_genes = set()
            for df in study_results:
                all_genes.update(df['gene'].tolist())
            
            # Create matrix of effect sizes and p-values
            effect_sizes = []
            pvalues = []
            studies = []
            
            for df in study_results:
                study_name = df['study'].iloc[0]
                studies.append(study_name)
                
                # Create series for this study
                study_effects = pd.Series(index=list(all_genes), dtype=float)
                study_pvals = pd.Series(index=list(all_genes), dtype=float)
                
                # Fill in values for genes present in this study
                for _, row in df.iterrows():
                    if row['gene'] in all_genes:
                        study_effects[row['gene']] = row['logfoldchanges']
                        study_pvals[row['gene']] = row['pvals_adj']
                
                effect_sizes.append(study_effects)
                pvalues.append(study_pvals)
            
            # Convert to DataFrame
            effects_df = pd.DataFrame(effect_sizes, index=studies).T
            pvals_df = pd.DataFrame(pvalues, index=studies).T
            
            # Remove genes with too many missing values
            min_studies = max(2, len(studies) // 2)
            valid_genes = effects_df.count(axis=1) >= min_studies
            effects_df = effects_df[valid_genes]
            pvals_df = pvals_df[valid_genes]
            
            # Perform meta-analysis using inverse variance weighting
            meta_results = []
            
            for gene in effects_df.index:
                gene_effects = effects_df.loc[gene].dropna()
                gene_pvals = pvals_df.loc[gene].dropna()
                
                if len(gene_effects) >= 2:
                    # Simple unweighted average (could be improved with proper meta-analysis)
                    meta_effect = gene_effects.mean()
                    meta_se = gene_effects.std() / np.sqrt(len(gene_effects))
                    
                    # Combine p-values using Fisher's method
                    combined_pval = stats.combine_pvalues(gene_pvals, method='fisher')[1]
                    
                    meta_results.append({
                        'gene': gene,
                        'meta_logfc': meta_effect,
                        'meta_se': meta_se,
                        'meta_pval': combined_pval,
                        'n_studies': len(gene_effects),
                        'heterogeneity': gene_effects.std()
                    })
            
            if meta_results:
                meta_df = pd.DataFrame(meta_results)
                
                # Multiple testing correction
                from scipy.stats import false_discovery_control
                meta_df['meta_pval_adj'] = false_discovery_control(meta_df['meta_pval'])
                
                # Sort by p-value
                meta_df = meta_df.sort_values('meta_pval')
                
                return meta_df
            
        except Exception as e:
            self.logger.error(f"Meta-analysis failed for {comparison_type}: {e}")
            
        return None
    
    def _meta_analyze_composition(self, datasets):
        """Meta-analyze cell type compositions across studies"""
        
        composition_results = {}
        
        # Calculate cell type proportions for each dataset
        proportions_list = []
        dataset_names = []
        
        for dataset_name, adata in datasets.items():
            if 'cell_type' in adata.obs.columns:
                proportions = adata.obs['cell_type'].value_counts(normalize=True)
                proportions_list.append(proportions)
                dataset_names.append(dataset_name)
        
        if len(proportions_list) > 1:
            # Combine proportions across datasets
            all_proportions = pd.DataFrame(proportions_list, index=dataset_names).fillna(0)
            
            # Calculate statistics
            composition_results['proportions'] = all_proportions
            composition_results['mean_proportions'] = all_proportions.mean()
            composition_results['std_proportions'] = all_proportions.std()
            composition_results['cv_proportions'] = all_proportions.std() / all_proportions.mean()
        
        return composition_results
    
    def _meta_analyze_expression(self, datasets):
        """Meta-analyze gene expression patterns"""
        
        expression_results = {}
        
        try:
            # Get common genes across all datasets
            common_genes = set()
            first = True
            
            for adata in datasets.values():
                if first:
                    common_genes = set(adata.var_names)
                    first = False
                else:
                    common_genes = common_genes.intersection(set(adata.var_names))
            
            if len(common_genes) > 100:  # Need sufficient overlap
                common_genes = list(common_genes)[:1000]  # Limit for computational efficiency
                
                # Calculate mean expression for each dataset
                dataset_expressions = []
                dataset_names = []
                
                for name, adata in datasets.items():
                    if 'cell_type' in adata.obs.columns:
                        # Calculate mean expression by cell type
                        cell_types = adata.obs['cell_type'].unique()
                        
                        for cell_type in cell_types:
                            cell_subset = adata[adata.obs['cell_type'] == cell_type]
                            if cell_subset.shape[0] > 10:  # Minimum cells
                                subset_genes = [g for g in common_genes if g in cell_subset.var_names]
                                if len(subset_genes) > 0:
                                    mean_exp = np.array(cell_subset[:, subset_genes].X.mean(axis=0)).flatten()
                                    dataset_expressions.append(mean_exp)
                                    dataset_names.append(f"{name}_{cell_type}")
                
                if len(dataset_expressions) > 1:
                    # Create expression matrix
                    expr_matrix = np.array(dataset_expressions)
                    
                    # Calculate correlations between datasets
                    correlation_matrix = np.corrcoef(expr_matrix)
                    
                    expression_results['correlation_matrix'] = correlation_matrix
                    expression_results['dataset_names'] = dataset_names
                    expression_results['common_genes_count'] = len(common_genes)
        
        except Exception as e:
            self.logger.warning(f"Expression meta-analysis failed: {e}")
        
        return expression_results
    
    def _calculate_heterogeneity(self, datasets):
        """Calculate heterogeneity metrics across studies"""
        
        heterogeneity = {}
        
        # Sample size heterogeneity
        sample_sizes = [adata.shape[0] for adata in datasets.values()]
        heterogeneity['sample_size_cv'] = np.std(sample_sizes) / np.mean(sample_sizes)
        
        # Gene count heterogeneity  
        gene_counts = [adata.shape[1] for adata in datasets.values()]
        heterogeneity['gene_count_cv'] = np.std(gene_counts) / np.mean(gene_counts)
        
        # Cell type diversity
        if all('cell_type' in adata.obs.columns for adata in datasets.values()):
            cell_type_counts = [len(adata.obs['cell_type'].unique()) for adata in datasets.values()]
            heterogeneity['cell_type_diversity_cv'] = np.std(cell_type_counts) / np.mean(cell_type_counts)
        
        return heterogeneity
    
    def _save_meta_results(self, results):
        """Save meta-analysis results"""
        
        output_dir = Path('results/tables')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        for analysis_name, analysis_results in results.items():
            if isinstance(analysis_results, pd.DataFrame):
                analysis_results.to_csv(output_dir / f'meta_{analysis_name}.csv', index=False)
            elif isinstance(analysis_results, dict):
                for key, value in analysis_results.items():
                    if isinstance(value, pd.DataFrame):
                        value.to_csv(output_dir / f'meta_{analysis_name}_{key}.csv')
                    elif isinstance(value, (pd.Series, np.ndarray)):
                        pd.DataFrame(value).to_csv(output_dir / f'meta_{analysis_name}_{key}.csv')
        
        self.logger.info("Meta-analysis results saved")
    
    def _plot_meta_results(self, results, datasets):
        """Create meta-analysis visualizations"""
        
        output_dir = Path('results/figures')
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Plot study overview
        self._plot_study_overview(datasets)
        
        # Plot composition meta-analysis
        if 'composition_meta' in results and results['composition_meta']:
            self._plot_composition_meta(results['composition_meta'])
        
        # Plot DE meta-analysis
        if 'de_meta' in results and results['de_meta']:
            self._plot_de_meta(results['de_meta'])
        
        # Plot expression correlations
        if 'expression_meta' in results and 'correlation_matrix' in results['expression_meta']:
            self._plot_expression_correlations(results['expression_meta'])
    
    def _plot_study_overview(self, datasets):
        """Plot overview of studies included in meta-analysis"""
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Sample sizes
        study_names = list(datasets.keys())
        sample_sizes = [adata.shape[0] for adata in datasets.values()]
        
        axes[0, 0].bar(study_names, sample_sizes)
        axes[0, 0].set_ylabel('Number of cells')
        axes[0, 0].set_title('Sample sizes across studies')
        axes[0, 0].tick_params(axis='x', rotation=45)
        
        # Gene counts
        gene_counts = [adata.shape[1] for adata in datasets.values()]
        axes[0, 1].bar(study_names, gene_counts)
        axes[0, 1].set_ylabel('Number of genes')
        axes[0, 1].set_title('Gene counts across studies')
        axes[0, 1].tick_params(axis='x', rotation=45)
        
        # Cell type diversity
        if all('cell_type' in adata.obs.columns for adata in datasets.values()):
            cell_type_counts = [len(adata.obs['cell_type'].unique()) for adata in datasets.values()]
            axes[1, 0].bar(study_names, cell_type_counts)
            axes[1, 0].set_ylabel('Number of cell types')
            axes[1, 0].set_title('Cell type diversity across studies')
            axes[1, 0].tick_params(axis='x', rotation=45)
        
        # Total UMI distribution
        total_umis = []
        study_labels = []
        for name, adata in datasets.items():
            umis = np.array(adata.X.sum(axis=1)).flatten()
            total_umis.extend(umis)
            study_labels.extend([name] * len(umis))
        
        # Box plot of UMI distributions
        import seaborn as sns
        plot_df = pd.DataFrame({'Study': study_labels, 'Total_UMI': total_umis})
        sns.boxplot(data=plot_df, x='Study', y='Total_UMI', ax=axes[1, 1])
        axes[1, 1].set_ylabel('Total UMI per cell')
        axes[1, 1].set_title('UMI distributions across studies')
        axes[1, 1].tick_params(axis='x', rotation=45)
        axes[1, 1].set_yscale('log')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'study_overview.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_composition_meta(self, composition_results):
        """Plot cell type composition meta-analysis"""
        
        if 'proportions' in composition_results:
            proportions_df = composition_results['proportions']
            
            fig, axes = plt.subplots(1, 2, figsize=(16, 6))
            
            # Heatmap of proportions
            sns.heatmap(proportions_df, annot=True, cmap='Blues', ax=axes[0])
            axes[0].set_title('Cell type proportions across studies')
            axes[0].set_ylabel('Study')
            axes[0].set_xlabel('Cell type')
            
            # Box plot of proportions
            proportions_melted = proportions_df.reset_index().melt(id_vars='index', 
                                                                   var_name='cell_type',
                                                                   value_name='proportion')
            proportions_melted.rename(columns={'index': 'study'}, inplace=True)
            
            sns.boxplot(data=proportions_melted, x='cell_type', y='proportion', ax=axes[1])
            axes[1].set_title('Cell type proportion distributions')
            axes[1].set_ylabel('Proportion')
            axes[1].tick_params(axis='x', rotation=45)
            
            plt.tight_layout()
            plt.savefig(output_dir / 'composition_meta_analysis.png', dpi=300, bbox_inches='tight')
            plt.close()
    
    def _plot_de_meta(self, de_meta_results):
        """Plot differential expression meta-analysis results"""
        
        for comparison, meta_df in list(de_meta_results.items())[:2]:  # Limit plots
            if isinstance(meta_df, pd.DataFrame) and len(meta_df) > 0:
                
                fig, axes = plt.subplots(1, 2, figsize=(16, 6))
                
                # Volcano plot
                meta_df['-log10_pval'] = -np.log10(meta_df['meta_pval_adj'].clip(lower=1e-300))
                
                colors = np.where(
                    (meta_df['meta_pval_adj'] < 0.05) & (np.abs(meta_df['meta_logfc']) > 0.5),
                    'red', 'gray'
                )
                
                axes[0].scatter(meta_df['meta_logfc'], meta_df['-log10_pval'], 
                               c=colors, alpha=0.6, s=10)
                axes[0].set_xlabel('Meta Log2 Fold Change')
                axes[0].set_ylabel('-Log10 Adjusted P-value')
                axes[0].set_title(f'Meta-analysis Volcano Plot: {comparison}')
                axes[0].axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
                axes[0].axvline(x=0.5, color='black', linestyle='--', alpha=0.5)
                axes[0].axvline(x=-0.5, color='black', linestyle='--', alpha=0.5)
                
                # Heterogeneity plot
                axes[1].scatter(meta_df['meta_logfc'], meta_df['heterogeneity'], alpha=0.6, s=10)
                axes[1].set_xlabel('Meta Log2 Fold Change')
                axes[1].set_ylabel('Heterogeneity (SD across studies)')
                axes[1].set_title(f'Effect Size Heterogeneity: {comparison}')
                
                plt.tight_layout()
                plt.savefig(output_dir / f'de_meta_analysis_{comparison}.png', 
                           dpi=300, bbox_inches='tight')
                plt.close()
    
    def _plot_expression_correlations(self, expression_results):
        """Plot expression correlation matrix"""
        
        correlation_matrix = expression_results['correlation_matrix']
        dataset_names = expression_results['dataset_names']
        
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Create correlation heatmap
        sns.heatmap(correlation_matrix, 
                    xticklabels=dataset_names,
                    yticklabels=dataset_names,
                    cmap='coolwarm', 
                    center=0,
                    annot=False,
                    ax=ax)
        
        ax.set_title('Expression Correlation Matrix Across Studies')
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        plt.tight_layout()
        plt.savefig(output_dir / 'expression_correlations.png', dpi=300, bbox_inches='tight')
        plt.close()
