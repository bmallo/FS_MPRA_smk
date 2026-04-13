#!/usr/bin/env python3
"""
Setup script for Fiberseq MPRA Analysis package
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="fiberseq_mpra",
    version="0.1.0",
    author="Fiberseq MPRA Analysis Team",
    description="Analysis pipeline for Fiber-seq MPRA footprint data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/your-org/fiberseq-mpra",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "pysam>=0.19.0",
        "pyyaml>=6.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=22.0.0",
            "flake8>=5.0.0",
        ],
        "viz": [
            "matplotlib>=3.5.0",
            "seaborn>=0.12.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "fiberseq-mpra=fiberseq_mpra.cli.main:main",
        ],
    },
)
