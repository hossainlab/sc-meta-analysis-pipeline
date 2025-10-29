"""
Test configuration and utilities
"""

import pytest
import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

@pytest.fixture
def test_config():
    """Provide test configuration"""
    return {
        'data': {
            'sources': ['manual'],
            'manual_files': [],
            'cache_dir': 'test_cache'
        },
        'quality_control': {
            'min_genes_per_cell': 100,
            'max_genes_per_cell': 3000,
            'min_cells_per_gene': 3,
            'max_mito_percent': 25
        },
        'logging': {
            'level': 'DEBUG',
            'file_logging': False,
            'console_logging': True
        }
    }

@pytest.fixture
def sample_adata():
    """Create sample AnnData object for testing"""
    import numpy as np
    import pandas as pd
    import anndata as ad
    
    # Create synthetic data
    n_cells = 100
    n_genes = 500
    
    X = np.random.poisson(2, size=(n_cells, n_genes))
    
    obs = pd.DataFrame({
        'cell_type': np.random.choice(['T cell', 'B cell', 'Monocyte'], n_cells),
        'condition': np.random.choice(['Control', 'Treatment'], n_cells)
    }, index=[f'Cell_{i}' for i in range(n_cells)])
    
    var = pd.DataFrame(
        index=[f'Gene_{i}' for i in range(n_genes)]
    )
    
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    return adata
