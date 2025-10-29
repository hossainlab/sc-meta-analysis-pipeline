"""
Setup script for Single Cell Meta-Analysis Pipeline
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="single-cell-meta-analysis",
    version="1.0.0",
    author="Single Cell Analysis Team",
    author_email="contact@example.com",
    description="A comprehensive pipeline for single-cell RNA-seq meta-analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/example/single-cell-meta-analysis",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
            "mypy>=0.812",
        ],
        "notebooks": [
            "jupyter>=1.0.0",
            "ipywidgets>=7.6.0",
            "papermill>=2.3.0",
        ]
    },
    entry_points={
        "console_scripts": [
            "sc-meta-analysis=run_pipeline:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
