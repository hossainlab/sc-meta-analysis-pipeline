"""
Data Collection Module
Handles downloading and loading single-cell datasets from various sources
"""

import os
import logging
import pandas as pd
import scanpy as sc
import anndata as ad
from pathlib import Path
import requests
from typing import Dict, List, Optional
import cellxgene_census

class DataCollector:
    """Handles data collection from various single-cell data sources"""
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.data_dir = Path(config.get('data_dir', 'data/raw'))
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
    def collect_datasets(self) -> Dict[str, ad.AnnData]:
        """Main method to collect all specified datasets"""
        datasets = {}
        
        # Collect from different sources
        if 'cellxgene' in self.config:
            datasets.update(self._collect_from_cellxgene())
        
        if 'geo' in self.config:
            datasets.update(self._collect_from_geo())
            
        if 'local_files' in self.config:
            datasets.update(self._collect_from_local())
            
        if 'urls' in self.config:
            datasets.update(self._collect_from_urls())
            
        return datasets
    
    def _collect_from_cellxgene(self) -> Dict[str, ad.AnnData]:
        """Collect datasets from CellxGene Census"""
        datasets = {}
        
        try:
            import cellxgene_census
            
            for dataset_info in self.config['cellxgene']:
                dataset_id = dataset_info['id']
                name = dataset_info.get('name', dataset_id)
                
                self.logger.info(f"Downloading dataset {name} from CellxGene...")
                
                # Download dataset
                adata = cellxgene_census.download_source_h5ad(dataset_id)
                
                # Add metadata
                adata.obs['study'] = name
                adata.obs['source'] = 'cellxgene'
                
                datasets[name] = adata
                
                # Save locally
                local_path = self.data_dir / f"{name}.h5ad"
                adata.write_h5ad(local_path)
                
                self.logger.info(f"Downloaded {name}: {adata.shape[0]} cells, {adata.shape[1]} genes")
                
        except Exception as e:
            self.logger.error(f"Error collecting from CellxGene: {str(e)}")
            
        return datasets
    
    def _collect_from_geo(self) -> Dict[str, ad.AnnData]:
        """Collect datasets from GEO"""
        datasets = {}
        
        try:
            from biomni.tool.database import query_geo
            
            for geo_info in self.config['geo']:
                accession = geo_info['accession']
                name = geo_info.get('name', accession)
                
                self.logger.info(f"Querying GEO for {accession}...")
                
                # Query GEO database
                geo_data = query_geo(accession)
                
                if geo_data and 'data not available' not in str(geo_data).lower():
                    # Process GEO data (this would need to be adapted based on actual GEO response)
                    self.logger.info(f"Found data for {accession}")
                    # Note: Actual implementation would depend on GEO data format
                else:
                    self.logger.warning(f"No data found for {accession}")
                    
        except Exception as e:
            self.logger.error(f"Error collecting from GEO: {str(e)}")
            
        return datasets
    
    def _collect_from_local(self) -> Dict[str, ad.AnnData]:
        """Load datasets from local files"""
        datasets = {}
        
        for file_info in self.config['local_files']:
            filepath = file_info['path']
            name = file_info.get('name', Path(filepath).stem)
            
            try:
                self.logger.info(f"Loading local file: {filepath}")
                
                if filepath.endswith('.h5ad'):
                    adata = sc.read_h5ad(filepath)
                elif filepath.endswith('.csv'):
                    df = pd.read_csv(filepath, index_col=0)
                    adata = ad.AnnData(df.T)  # Transpose so cells are rows
                elif filepath.endswith('.xlsx'):
                    df = pd.read_excel(filepath, index_col=0)
                    adata = ad.AnnData(df.T)
                else:
                    self.logger.warning(f"Unsupported file format: {filepath}")
                    continue
                
                # Add metadata
                adata.obs['study'] = name
                adata.obs['source'] = 'local'
                
                datasets[name] = adata
                self.logger.info(f"Loaded {name}: {adata.shape[0]} cells, {adata.shape[1]} genes")
                
            except Exception as e:
                self.logger.error(f"Error loading {filepath}: {str(e)}")
                
        return datasets
    
    def _collect_from_urls(self) -> Dict[str, ad.AnnData]:
        """Download datasets from URLs"""
        datasets = {}
        
        for url_info in self.config['urls']:
            url = url_info['url']
            name = url_info.get('name', Path(url).stem)
            
            try:
                self.logger.info(f"Downloading from URL: {url}")
                
                # Download file
                local_path = self.data_dir / Path(url).name
                self._download_file(url, local_path)
                
                # Load based on file extension
                if str(local_path).endswith('.h5ad'):
                    adata = sc.read_h5ad(local_path)
                elif str(local_path).endswith('.csv'):
                    df = pd.read_csv(local_path, index_col=0)
                    adata = ad.AnnData(df.T)
                else:
                    self.logger.warning(f"Unsupported file format: {local_path}")
                    continue
                
                # Add metadata
                adata.obs['study'] = name
                adata.obs['source'] = 'url'
                
                datasets[name] = adata
                self.logger.info(f"Downloaded {name}: {adata.shape[0]} cells, {adata.shape[1]} genes")
                
            except Exception as e:
                self.logger.error(f"Error downloading {url}: {str(e)}")
                
        return datasets
    
    def _download_file(self, url: str, local_path: Path):
        """Download file from URL"""
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        with open(local_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
    
    def validate_datasets(self, datasets: Dict[str, ad.AnnData]) -> Dict[str, ad.AnnData]:
        """Validate and clean datasets"""
        validated_datasets = {}
        
        for name, adata in datasets.items():
            try:
                # Basic validation
                if adata.shape[0] == 0 or adata.shape[1] == 0:
                    self.logger.warning(f"Empty dataset: {name}")
                    continue
                
                # Make variable names unique
                adata.var_names_make_unique()
                
                # Ensure X is not None
                if adata.X is None:
                    self.logger.warning(f"No expression data in {name}")
                    continue
                
                validated_datasets[name] = adata
                
            except Exception as e:
                self.logger.error(f"Validation failed for {name}: {str(e)}")
                
        return validated_datasets
