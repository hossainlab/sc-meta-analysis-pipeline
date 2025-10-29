# Troubleshooting Guide

## Common Issues and Solutions

### Installation Issues

**Issue**: Package installation fails
```bash
ERROR: Could not install packages due to an EnvironmentError
```

**Solutions**:
1. Create a fresh conda environment:
```bash
conda create -n scmeta python=3.9
conda activate scmeta
pip install -r requirements.txt
```

2. Use mamba for faster dependency resolution:
```bash
conda install mamba
mamba install -c conda-forge -c bioconda scanpy scvi-tools
```

### Data Loading Issues

**Issue**: "File not found" error
```
FileNotFoundError: File not found: data/raw/dataset.h5ad
```

**Solutions**:
1. Check file paths in configuration
2. Ensure data files are in correct directory structure
3. Use absolute paths if necessary

**Issue**: Memory errors when loading large datasets
```
MemoryError: Unable to allocate array
```

**Solutions**:
1. Increase memory limits in config
2. Use chunked processing for large files
3. Consider using backed mode: `adata = sc.read_h5ad('file.h5ad', backed='r')`

### Quality Control Issues

**Issue**: Too few cells remaining after filtering
```
WARNING: Only 50 cells remaining after QC filtering
```

**Solutions**:
1. Relax filtering parameters in config:
```yaml
quality_control:
  min_genes_per_cell: 100  # Reduce from 200
  max_mito_percent: 30     # Increase from 20
```

2. Check data quality before filtering
3. Adjust thresholds based on tissue type

### Normalization Issues

**Issue**: HVG selection fails
```
ValueError: No highly variable genes found
```

**Solutions**:
1. Check if data is already log-transformed
2. Try different HVG methods:
```yaml
normalization:
  hvg_method: seurat  # Try 'seurat' instead of 'seurat_v3'
  n_top_genes: 1000   # Reduce number
```

### Integration Issues

**Issue**: Harmony integration fails
```
ModuleNotFoundError: No module named 'harmony'
```

**Solutions**:
1. Install harmony:
```bash
pip install harmonypy
```

2. Try alternative methods:
```yaml
batch_correction:
  method: combat  # or scanorama
```

**Issue**: Poor integration quality
```
WARNING: High silhouette score for batches indicates poor mixing
```

**Solutions**:
1. Try different integration methods
2. Adjust parameters:
```yaml
batch_correction:
  method: scvi  # More sophisticated method
```

3. Check if batch correction is necessary

### Cell Type Annotation Issues

**Issue**: Annotation confidence too low
```
WARNING: Low confidence cell type annotations
```

**Solutions**:
1. Lower confidence threshold:
```yaml
annotation:
  confidence_threshold: 0.3  # Reduce from 0.5
```

2. Use manual annotation with known markers
3. Try reference-based annotation

### Differential Expression Issues

**Issue**: No significant genes found
```
WARNING: No genes pass significance threshold
```

**Solutions**:
1. Relax thresholds:
```yaml
differential_expression:
  min_pct: 0.05        # Reduce from 0.1
  logfc_threshold: 0.1  # Reduce from 0.25
```

2. Check sample sizes per group
3. Try different statistical methods

### Visualization Issues

**Issue**: Plots are empty or corrupted
```
WARNING: Could not create UMAP plot
```

**Solutions**:
1. Ensure UMAP coordinates exist:
```python
sc.tl.umap(adata)  # Compute UMAP if missing
```

2. Check figure directory permissions
3. Try different plot formats

### Performance Issues

**Issue**: Pipeline runs very slowly
```
INFO: Processing taking longer than expected
```

**Solutions**:
1. Increase CPU cores:
```yaml
performance:
  n_cores: 8  # Use more cores
```

2. Use GPU acceleration if available:
```yaml
performance:
  use_gpu: true
```

3. Process datasets separately then integrate

**Issue**: Out of memory errors
```
MemoryError: Unable to allocate array
```

**Solutions**:
1. Increase memory limit:
```yaml
performance:
  memory_limit_gb: 32  # Increase from 16
```

2. Use chunked processing
3. Filter data more aggressively before analysis

## Debug Mode

Run pipeline in debug mode for detailed information:

```bash
python run_pipeline.py --config config/debug_config.yaml --debug
```

Create debug config with verbose logging:
```yaml
logging:
  level: DEBUG
  file_logging: true
  console_logging: true
```

## Getting Help

1. Check log files in `logs/pipeline.log`
2. Enable debug mode for detailed output
3. Review configuration settings
4. Check data integrity before processing
5. Try with smaller test datasets first

## Common Configuration Fixes

### For PBMC datasets:
```yaml
quality_control:
  max_mito_percent: 20
  
normalization:
  hvg_method: seurat_v3
  n_top_genes: 2000
```

### For tissue-specific datasets:
```yaml
quality_control:
  max_mito_percent: 30  # Higher for some tissues
  
batch_correction:
  method: harmony  # Good for tissue integration
```

### For large-scale studies:
```yaml
performance:
  n_cores: 16
  memory_limit_gb: 64
  chunk_size: 2000
  
differential_expression:
  max_cells_per_group: 500  # Downsample for speed
```
