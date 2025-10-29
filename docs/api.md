# API Documentation

## Core Modules

### data_collection.DataCollector

Main class for collecting and loading single-cell datasets.

```python
from data_collection import DataCollector

config = {'sources': ['geo'], 'geo_datasets': ['GSE123456']}
collector = DataCollector(config)
datasets = collector.collect_datasets()
```

**Methods:**
- `collect_datasets()`: Collect all configured datasets
- `load_geo_dataset(accession)`: Load specific GEO dataset
- `load_manual_files(file_paths)`: Load manual data files

### quality_control.QualityController

Handles quality control and filtering of single-cell data.

```python
from quality_control import QualityController

qc = QualityController(config['quality_control'])
adata_filtered = qc.filter_cells_and_genes(adata)
```

**Methods:**
- `filter_cells_and_genes(adata)`: Apply cell and gene filtering
- `detect_doublets(adata)`: Detect and optionally remove doublets
- `calculate_qc_metrics(adata)`: Calculate quality control metrics

### normalization.Normalizer

Implements normalization and feature selection.

```python
from normalization import Normalizer

normalizer = Normalizer(config['normalization'])
adata_norm = normalizer.normalize(adata)
```

**Methods:**
- `normalize(adata)`: Apply normalization pipeline
- `find_highly_variable_genes(adata)`: Identify highly variable genes
- `calculate_size_factors(adata)`: Calculate normalization factors

### batch_correction.BatchCorrector

Handles batch correction and dataset integration.

```python
from batch_correction import BatchCorrector

corrector = BatchCorrector(config['batch_correction'])
adata_integrated = corrector.integrate_datasets(datasets)
```

**Methods:**
- `integrate_datasets(datasets)`: Integrate multiple datasets
- `calculate_integration_metrics(adata)`: Assess integration quality

### differential_expression.DifferentialExpressionAnalyzer

Performs differential expression analysis.

```python
from differential_expression import DifferentialExpressionAnalyzer

de_analyzer = DifferentialExpressionAnalyzer(config['differential_expression'])
de_results = de_analyzer.find_markers(adata)
```

**Methods:**
- `find_markers(adata)`: Find cell type marker genes
- `find_condition_de_genes(adata)`: Find condition-specific DE genes

## Utility Functions

### utils.config

Configuration management utilities.

```python
from utils.config import load_config, setup_logging

config = load_config('config/my_config.yaml')
logger = setup_logging(config)
```

### utils.file_utils

File handling utilities.

```python
from utils.file_utils import discover_data_files, load_data_file

files = discover_data_files('data/raw')
adata = load_data_file('data/raw/dataset.h5ad')
```

### utils.visualization

Visualization utilities.

```python
from utils.visualization import create_qc_dashboard, create_integration_plot

create_qc_dashboard(adata, 'results/figures')
create_integration_plot(adata, 'results/figures')
```

## Configuration Schema

The pipeline uses YAML configuration files with the following structure:

```yaml
data:
  sources: [geo, manual]
  geo_datasets: []
  manual_files: []
  species: human

quality_control:
  min_genes_per_cell: 200
  max_genes_per_cell: 5000
  max_mito_percent: 20
  doublet_detection: true

normalization:
  method: total_count
  target_sum: 10000
  hvg_method: seurat_v3
  n_top_genes: 2000

batch_correction:
  method: harmony
  batch_key: dataset

annotation:
  method: automatic
  use_markers_db: true

differential_expression:
  method: wilcoxon
  min_pct: 0.1
  logfc_threshold: 0.25

enrichment:
  gene_sets:
    - MSigDB_Hallmark
    - MSigDB_C2_CP_REACTOME
  fdr_threshold: 0.05

visualization:
  create_qc_plots: true
  create_umap_plots: true
  dpi: 300
```
