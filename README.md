# Single Cell Meta-Analysis Pipeline

A comprehensive, reproducible pipeline for single-cell RNA-seq meta-analysis using modern computational tools and best practices.

## 🚀 Features

- **Data Collection**: Automated downloading from GEO, manual data import
- **Quality Control**: Cell/gene filtering, doublet detection, QC visualization
- **Normalization**: Multiple normalization methods (total count, scran, etc.)
- **Batch Correction**: Harmony, Scanorama, scVI, ComBat integration methods
- **Cell Annotation**: Automated and manual cell type annotation
- **Differential Expression**: Multiple DE methods with meta-analysis
- **Pathway Enrichment**: GSEA with MSigDB and custom gene sets
- **Meta-Analysis**: Statistical combination across multiple studies
- **Visualization**: Publication-ready plots and interactive visualizations
- **Reproducibility**: Version-controlled, configurable, fully documented

## 📋 Requirements

- Python 3.8+
- 8GB+ RAM recommended
- 50GB+ disk space for large datasets

## 🔧 Installation

### Option 1: Standard Installation

```bash
# Clone the repository
git clone https://github.com/example/single-cell-meta-analysis.git
cd single-cell-meta-analysis

# Create conda environment
conda create -n sc-meta python=3.9
conda activate sc-meta

# Install dependencies
pip install -r requirements.txt

# Install the package
pip install -e .
```

### Option 2: Docker Installation

```bash
# Build Docker image
docker build -t sc-meta-analysis .

# Run container
docker run -it -v $(pwd)/data:/app/data -v $(pwd)/results:/app/results sc-meta-analysis
```

## 🚦 Quick Start

### 1. Basic Usage

```bash
# Run complete pipeline with default configuration
python run_pipeline.py --config config/default_config.yaml

# Run specific analysis step
python run_pipeline.py --config config/default_config.yaml --step quality_control
```

### 2. Custom Configuration

Create a custom config file:

```yaml
# my_config.yaml
data:
  geo_datasets:
    - GSE123456
    - GSE789012
  manual_files:
    - data/raw/dataset1.h5ad
    - data/raw/dataset2.h5ad

quality_control:
  min_genes_per_cell: 300
  max_mito_percent: 15

batch_correction:
  method: harmony

annotation:
  method: automatic
```

Run with custom config:
```bash
python run_pipeline.py --config my_config.yaml
```

### 3. Python API Usage

```python
import sys
sys.path.insert(0, 'src')

from data_collection import DataCollector
from quality_control import QualityController
from utils.config import load_config

# Load configuration
config = load_config('config/default_config.yaml')

# Initialize components
collector = DataCollector(config['data'])
qc = QualityController(config['quality_control'])

# Collect data
datasets = collector.collect_datasets()

# Run quality control
filtered_datasets = {}
for name, adata in datasets.items():
    filtered_datasets[name] = qc.filter_cells_and_genes(adata)
```

## 📁 Pipeline Structure

```
single_cell_meta_analysis_pipeline/
├── run_pipeline.py          # Main pipeline script
├── setup.py                 # Package setup
├── requirements.txt         # Dependencies
├── README.md               # This file
├── LICENSE                 # License file
├── config/                 # Configuration files
│   └── default_config.yaml
├── src/                    # Source code
│   ├── data_collection.py
│   ├── quality_control.py
│   ├── normalization.py
│   ├── batch_correction.py
│   ├── cell_annotation.py
│   ├── differential_expression.py
│   ├── enrichment_analysis.py
│   ├── meta_analysis.py
│   ├── visualization.py
│   └── utils/
│       ├── config.py
│       ├── logger.py
│       └── data_utils.py
├── data/                   # Data directory
│   ├── raw/               # Raw input data
│   └── processed/         # Processed data
├── results/               # Analysis results
│   ├── figures/          # Generated plots
│   └── tables/           # Result tables
├── notebooks/            # Jupyter notebooks
│   ├── 01_data_exploration.ipynb
│   ├── 02_quality_control.ipynb
│   ├── 03_integration_analysis.ipynb
│   └── 04_meta_analysis.ipynb
├── tests/               # Unit tests
├── docs/               # Documentation
└── examples/           # Example datasets and configs
```

## 🔬 Analysis Steps

### Step 1: Data Collection
- Download datasets from GEO
- Load manual data files
- Convert formats to AnnData
- Basic metadata extraction

### Step 2: Quality Control
- Filter low-quality cells and genes
- Detect and remove doublets
- Calculate QC metrics
- Generate QC visualizations

### Step 3: Normalization
- Normalize total counts
- Log transformation
- Identify highly variable genes
- Scale data

### Step 4: Batch Correction
- Integrate multiple datasets
- Remove batch effects
- Preserve biological variation
- Generate integration plots

### Step 5: Cell Type Annotation
- Automated annotation using markers
- Manual annotation support
- Transfer labels between datasets
- Validate annotations

### Step 6: Differential Expression
- Find marker genes for cell types
- Compare conditions/treatments
- Support multiple DE methods
- Generate volcano plots

### Step 7: Pathway Enrichment
- GSEA analysis
- Multiple gene set databases
- Pathway visualization
- Functional interpretation

### Step 8: Meta-Analysis
- Combine results across studies
- Statistical meta-analysis
- Heterogeneity assessment
- Effect size estimation

### Step 9: Visualization
- UMAP/t-SNE plots
- Heatmaps and dot plots
- Volcano plots
- Publication-ready figures

## 📊 Supported Data Formats

### Input Formats
- **H5AD**: AnnData HDF5 format (recommended)
- **CSV/TSV**: Gene expression matrices
- **MTX**: Matrix Market format (10X Genomics)
- **H5**: HDF5 files
- **Excel**: XLSX/XLS files

### Output Formats
- **H5AD**: Processed AnnData objects
- **CSV**: Result tables
- **PNG/PDF**: High-resolution figures
- **HTML**: Interactive plots

## ⚙️ Configuration Options

### Data Collection
```yaml
data:
  sources: [geo, manual]
  geo_datasets: [GSE123456]
  manual_files: [data/dataset.h5ad]
  cache_dir: data/cache
  species: human
```

### Quality Control
```yaml
quality_control:
  min_genes_per_cell: 200
  max_genes_per_cell: 5000
  min_cells_per_gene: 3
  max_mito_percent: 20
  doublet_detection: true
```

### Batch Correction
```yaml
batch_correction:
  method: harmony  # harmony, scanorama, scvi, combat
  batch_key: dataset
  use_pca: true
  n_components: 50
```

## 📈 Output Files

### Tables
- `quality_control_metrics.csv`: QC statistics
- `de_cell_type_markers_*.csv`: Marker genes
- `enrichment_results_*.csv`: Pathway analysis
- `meta_analysis_results.csv`: Combined results

### Figures
- `quality_control_overview.png`: QC metrics
- `integration_results.png`: Batch correction
- `cell_type_umap.png`: Annotated UMAP
- `marker_genes_heatmap.png`: Expression heatmap
- `enrichment_dotplot.png`: Pathway results

## 🧪 Testing

Run unit tests:
```bash
pytest tests/
```

Run integration tests:
```bash
pytest tests/test_integration.py -v
```

## 🐛 Troubleshooting

### Common Issues

1. **Memory Errors**
   - Increase memory limit in config
   - Use data chunking
   - Filter datasets to smaller size

2. **Missing Dependencies**
   ```bash
   pip install -r requirements.txt
   conda install -c bioconda scanpy
   ```

3. **R Dependencies** (for scran normalization)
   ```bash
   conda install -c conda-forge r-base
   pip install rpy2
   ```

4. **GPU Support** (for scVI)
   ```bash
   pip install scvi-tools[cuda]
   ```

### Performance Optimization

- Enable GPU acceleration for scVI
- Increase number of cores in config
- Use SSD storage for data
- Pre-filter large datasets

## 📚 Documentation

- [Full API Documentation](docs/api.md)
- [Tutorial Notebooks](notebooks/)
- [Best Practices Guide](docs/best_practices.md)
- [Troubleshooting Guide](docs/troubleshooting.md)

## 🤝 Contributing

1. Fork the repository
2. Create feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open Pull Request

## 📄 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file.

## 📞 Support

- Create GitHub issue for bugs/features
- Email: support@example.com
- Documentation: https://docs.example.com

## 🔬 Citation

If you use this pipeline in your research, please cite:

```
Single Cell Meta-Analysis Pipeline (2024)
Available at: https://github.com/example/single-cell-meta-analysis
```

## 🙏 Acknowledgments

- scanpy developers for single-cell analysis tools
- scVI-tools team for probabilistic models
- Decoupler team for pathway analysis
- The single-cell genomics community

---

## 🗺️ Roadmap

### Version 1.1 (Coming Soon)
- [ ] Spatial transcriptomics support
- [ ] Multi-modal analysis (CITE-seq, ATAC-seq)
- [ ] Real-time analysis dashboard
- [ ] Cloud deployment options

### Version 1.2 (Future)
- [ ] Machine learning cell type prediction
- [ ] Drug target identification
- [ ] Trajectory analysis integration
- [ ] Automated report generation

---

**Happy analyzing! 🧬✨**
