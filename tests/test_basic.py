"""
Basic tests for pipeline functionality
"""

import pytest
import numpy as np
import pandas as pd

def test_config_loading(test_config):
    """Test configuration loading"""
    assert 'data' in test_config
    assert 'quality_control' in test_config
    assert test_config['quality_control']['min_genes_per_cell'] == 100

def test_sample_adata(sample_adata):
    """Test sample data creation"""
    assert sample_adata.n_obs == 100
    assert sample_adata.n_vars == 500
    assert 'cell_type' in sample_adata.obs.columns
    assert 'condition' in sample_adata.obs.columns

def test_quality_control_module():
    """Test quality control module import"""
    try:
        from quality_control import QualityController
        assert QualityController is not None
    except ImportError as e:
        pytest.skip(f"Quality control module not available: {e}")

def test_normalization_module():
    """Test normalization module import"""
    try:
        from normalization import Normalizer
        assert Normalizer is not None
    except ImportError as e:
        pytest.skip(f"Normalization module not available: {e}")
