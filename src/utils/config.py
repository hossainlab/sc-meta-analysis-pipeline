"""
Configuration utilities for the pipeline
"""

import yaml
import logging
from pathlib import Path
from typing import Dict, Any

def load_config(config_path: str) -> Dict[str, Any]:
    """Load configuration from YAML file"""
    
    config_file = Path(config_path)
    
    if not config_file.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    # Validate required sections
    required_sections = [
        'data', 'quality_control', 'normalization', 
        'batch_correction', 'annotation', 'differential_expression'
    ]
    
    for section in required_sections:
        if section not in config:
            raise ValueError(f"Required configuration section missing: {section}")
    
    return config

def setup_logging(config: Dict[str, Any]) -> logging.Logger:
    """Setup logging configuration"""
    
    logging_config = config.get('logging', {})
    level = logging_config.get('level', 'INFO')
    file_logging = logging_config.get('file_logging', True)
    console_logging = logging_config.get('console_logging', True)
    
    # Create logger
    logger = logging.getLogger('single_cell_pipeline')
    logger.setLevel(getattr(logging, level.upper()))
    
    # Clear existing handlers
    logger.handlers = []
    
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Console handler
    if console_logging:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    
    # File handler
    if file_logging:
        log_dir = Path('logs')
        log_dir.mkdir(exist_ok=True)
        
        file_handler = logging.FileHandler(log_dir / 'pipeline.log')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger

def validate_config(config: Dict[str, Any]) -> bool:
    """Validate configuration values"""
    
    # Validate QC parameters
    qc_config = config['quality_control']
    
    if qc_config['min_genes_per_cell'] >= qc_config['max_genes_per_cell']:
        raise ValueError("min_genes_per_cell must be less than max_genes_per_cell")
    
    if qc_config['max_mito_percent'] < 0 or qc_config['max_mito_percent'] > 100:
        raise ValueError("max_mito_percent must be between 0 and 100")
    
    # Validate normalization parameters
    norm_config = config['normalization']
    
    if norm_config['target_sum'] <= 0:
        raise ValueError("target_sum must be positive")
    
    if norm_config['n_top_genes'] <= 0:
        raise ValueError("n_top_genes must be positive")
    
    # Validate DE parameters
    de_config = config['differential_expression']
    
    if de_config['min_pct'] < 0 or de_config['min_pct'] > 1:
        raise ValueError("min_pct must be between 0 and 1")
    
    if de_config['logfc_threshold'] < 0:
        raise ValueError("logfc_threshold must be non-negative")
    
    return True

def create_output_dirs(config: Dict[str, Any]) -> Dict[str, Path]:
    """Create necessary output directories"""
    
    dirs = {
        'results': Path('results'),
        'figures': Path('results/figures'),
        'tables': Path('results/tables'),
        'processed_data': Path('data/processed'),
        'cache': Path(config['data'].get('cache_dir', 'data/cache')),
        'logs': Path('logs')
    }
    
    for name, path in dirs.items():
        path.mkdir(parents=True, exist_ok=True)
    
    return dirs

def get_available_methods() -> Dict[str, list]:
    """Get available methods for each analysis step"""
    
    return {
        'normalization': ['total_count', 'median', 'scran'],
        'batch_correction': ['harmony', 'scanorama', 'scvi', 'combat'],
        'annotation': ['automatic', 'manual', 'reference'],
        'differential_expression': ['wilcoxon', 't-test', 'logreg'],
        'enrichment': ['gsea', 'over_representation', 'aucell']
    }
