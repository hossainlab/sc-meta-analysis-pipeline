#!/usr/bin/env python3
"""
Single Cell Meta-Analysis Pipeline
Main execution script for the complete analysis workflow
"""

import os
import sys
import argparse
import logging
from pathlib import Path
import yaml

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from data_collection import DataCollector
from quality_control import QualityController
from normalization import Normalizer
from batch_correction import BatchCorrector
from cell_annotation import CellAnnotator
from differential_expression import DifferentialExpressionAnalyzer
from enrichment_analysis import EnrichmentAnalyzer
from meta_analysis import MetaAnalyzer
from visualization import Visualizer
from utils.logger import setup_logging
from utils.config import load_config

class SingleCellMetaAnalysisPipeline:
    """Main pipeline class for single-cell meta-analysis"""
    
    def __init__(self, config_path):
        self.config = load_config(config_path)
        self.setup_directories()
        setup_logging(self.config['logging'])
        self.logger = logging.getLogger(__name__)
        
    def setup_directories(self):
        """Create necessary output directories"""
        for dir_name in ['data/processed', 'results/figures', 'results/tables']:
            Path(dir_name).mkdir(parents=True, exist_ok=True)
    
    def run_data_collection(self):
        """Step 1: Data Collection"""
        self.logger.info("Starting data collection...")
        collector = DataCollector(self.config['data'])
        datasets = collector.collect_datasets()
        self.logger.info(f"Collected {len(datasets)} datasets")
        return datasets
    
    def run_quality_control(self, datasets):
        """Step 2: Quality Control"""
        self.logger.info("Starting quality control...")
        qc = QualityController(self.config['quality_control'])
        filtered_datasets = {}
        for name, adata in datasets.items():
            filtered_datasets[name] = qc.filter_cells_and_genes(adata)
        return filtered_datasets
    
    def run_normalization(self, datasets):
        """Step 3: Normalization"""
        self.logger.info("Starting normalization...")
        normalizer = Normalizer(self.config['normalization'])
        normalized_datasets = {}
        for name, adata in datasets.items():
            normalized_datasets[name] = normalizer.normalize(adata)
        return normalized_datasets
    
    def run_batch_correction(self, datasets):
        """Step 4: Batch Correction"""
        self.logger.info("Starting batch correction...")
        corrector = BatchCorrector(self.config['batch_correction'])
        corrected_data = corrector.integrate_datasets(datasets)
        return corrected_data
    
    def run_cell_annotation(self, adata):
        """Step 5: Cell Type Annotation"""
        self.logger.info("Starting cell annotation...")
        annotator = CellAnnotator(self.config['annotation'])
        annotated_data = annotator.annotate_cells(adata)
        return annotated_data
    
    def run_differential_expression(self, adata):
        """Step 6: Differential Expression Analysis"""
        self.logger.info("Starting differential expression analysis...")
        de_analyzer = DifferentialExpressionAnalyzer(self.config['differential_expression'])
        de_results = de_analyzer.find_markers(adata)
        return de_results
    
    def run_enrichment_analysis(self, de_results, adata):
        """Step 7: Pathway Enrichment Analysis"""
        self.logger.info("Starting enrichment analysis...")
        enrichment_analyzer = EnrichmentAnalyzer(self.config['enrichment'])
        enrichment_results = enrichment_analyzer.run_gsea(de_results, adata)
        return enrichment_results
    
    def run_meta_analysis(self, datasets, de_results):
        """Step 8: Meta-Analysis"""
        self.logger.info("Starting meta-analysis...")
        meta_analyzer = MetaAnalyzer(self.config['meta_analysis'])
        meta_results = meta_analyzer.combine_studies(datasets, de_results)
        return meta_results
    
    def run_visualization(self, adata, de_results, enrichment_results, meta_results):
        """Step 9: Visualization"""
        self.logger.info("Creating visualizations...")
        visualizer = Visualizer(self.config['visualization'])
        visualizer.create_all_plots(adata, de_results, enrichment_results, meta_results)
    
    def run_complete_pipeline(self):
        """Run the complete analysis pipeline"""
        try:
            # Step 1: Data Collection
            datasets = self.run_data_collection()
            
            # Step 2: Quality Control
            datasets = self.run_quality_control(datasets)
            
            # Step 3: Normalization
            datasets = self.run_normalization(datasets)
            
            # Step 4: Batch Correction
            integrated_data = self.run_batch_correction(datasets)
            
            # Step 5: Cell Annotation
            annotated_data = self.run_cell_annotation(integrated_data)
            
            # Step 6: Differential Expression
            de_results = self.run_differential_expression(annotated_data)
            
            # Step 7: Enrichment Analysis
            enrichment_results = self.run_enrichment_analysis(de_results, annotated_data)
            
            # Step 8: Meta-Analysis
            meta_results = self.run_meta_analysis(datasets, de_results)
            
            # Step 9: Visualization
            self.run_visualization(annotated_data, de_results, enrichment_results, meta_results)
            
            self.logger.info("Pipeline completed successfully!")
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            raise

def main():
    parser = argparse.ArgumentParser(description='Single Cell Meta-Analysis Pipeline')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--step', help='Run specific step only', 
                       choices=['data_collection', 'quality_control', 'normalization', 
                               'batch_correction', 'annotation', 'differential_expression',
                               'enrichment', 'meta_analysis', 'visualization'])
    
    args = parser.parse_args()
    
    pipeline = SingleCellMetaAnalysisPipeline(args.config)
    
    if args.step:
        # Run specific step
        step_methods = {
            'data_collection': pipeline.run_data_collection,
            'quality_control': lambda: pipeline.run_quality_control({}),
            'normalization': lambda: pipeline.run_normalization({}),
            'batch_correction': lambda: pipeline.run_batch_correction({}),
            'annotation': lambda: pipeline.run_cell_annotation(None),
            'differential_expression': lambda: pipeline.run_differential_expression(None),
            'enrichment': lambda: pipeline.run_enrichment_analysis({}, None),
            'meta_analysis': lambda: pipeline.run_meta_analysis({}, {}),
            'visualization': lambda: pipeline.run_visualization(None, {}, {}, {})
        }
        step_methods[args.step]()
    else:
        # Run complete pipeline
        pipeline.run_complete_pipeline()

if __name__ == "__main__":
    main()
