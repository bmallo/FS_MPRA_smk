#!/usr/bin/env python3
"""
Configuration Management for Fiberseq MPRA Analysis

This module handles configuration loading, validation, and defaults.
Configuration can be provided via YAML files, command-line arguments, or programmatically.
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import yaml

logger = logging.getLogger(__name__)


# Default configuration values
DEFAULT_CONFIG = {
    # Filtering parameters
    'min_variant_reads': 500,
    'nucleosome_range': None,  # None means no filtering, or [min, max]
    'require_single_variant': True,
    'min_read_length': 0,
    
    # Footprint size bins
    'size_bins': {
        'TF': [20, 49],
        'mid': [50, 79],
        'nucleosome': [80, 200],
    },
    
    # Analysis region
    'position_range': None,  # None means use full read length, or [start, end]
    
    # Statistical parameters
    'fdr_threshold': 0.05,
    'min_effect_size': 0.0,  # Minimum |log2FC| to report
    
    # Output parameters
    'output_format': ['html', 'tsv'],
    'per_variant_output': False,
    'generate_static_figures': True,
    
    # Visualization parameters
    'heatmap_cmap': 'RdBu_r',
    'heatmap_vmin': -3,  # log2FC scale
    'heatmap_vmax': 3,
    'significance_cmap': 'viridis',
    'significance_vmax': 10,  # -log10(pvalue) scale
}


@dataclass
class Config:
    """
    Configuration container for Fiberseq MPRA analysis.
    
    Attributes are organized by category:
    - Filtering: Read and variant filtering parameters
    - Size bins: Footprint size bin definitions
    - Analysis: Region and statistical parameters
    - Output: Output format and path settings
    - Visualization: Plot styling parameters
    """
    
    # Filtering parameters
    min_variant_reads: int = 500
    nucleosome_range: Optional[Tuple[int, int]] = None
    require_single_variant: bool = True
    min_read_length: int = 0
    
    # Footprint size bins (name -> [min, max])
    size_bins: Dict[str, List[int]] = field(default_factory=lambda: {
        'TF': [20, 49],
        'mid': [50, 79],
        'nucleosome': [80, 200],
    })
    
    # Analysis region
    position_range: Optional[Tuple[int, int]] = None
    
    # Statistical parameters
    fdr_threshold: float = 0.05
    min_effect_size: float = 0.0
    
    # Output parameters
    output_format: List[str] = field(default_factory=lambda: ['html', 'tsv'])
    per_variant_output: bool = False
    generate_static_figures: bool = True
    
    # Visualization parameters
    heatmap_cmap: str = 'RdBu_r'
    heatmap_vmin: float = -3.0
    heatmap_vmax: float = 3.0
    significance_cmap: str = 'viridis'
    significance_vmax: float = 10.0
    
    def validate(self) -> List[str]:
        """
        Validate configuration values.
        
        Returns:
        --------
        List[str]
            List of validation error messages (empty if valid)
        """
        errors = []
        
        # Check min_variant_reads
        if self.min_variant_reads < 1:
            errors.append("min_variant_reads must be >= 1")
        
        # Check nucleosome_range
        if self.nucleosome_range is not None:
            if len(self.nucleosome_range) != 2:
                errors.append("nucleosome_range must be [min, max] or None")
            elif self.nucleosome_range[0] > self.nucleosome_range[1]:
                errors.append("nucleosome_range min must be <= max")
        
        # Check size_bins
        if not self.size_bins:
            errors.append("size_bins cannot be empty")
        for name, (min_size, max_size) in self.size_bins.items():
            if min_size > max_size:
                errors.append(f"size_bin '{name}': min ({min_size}) > max ({max_size})")
            if min_size < 1:
                errors.append(f"size_bin '{name}': min must be >= 1")
        
        # Check position_range
        if self.position_range is not None:
            if len(self.position_range) != 2:
                errors.append("position_range must be [start, end] or None")
            elif self.position_range[0] > self.position_range[1]:
                errors.append("position_range start must be <= end")
        
        # Check fdr_threshold
        if not 0 < self.fdr_threshold < 1:
            errors.append("fdr_threshold must be between 0 and 1")
        
        # Check output_format
        valid_formats = {'html', 'tsv', 'pdf', 'json'}
        for fmt in self.output_format:
            if fmt not in valid_formats:
                errors.append(f"Unknown output format: {fmt}")
        
        return errors
    
    def to_dict(self) -> Dict:
        """Convert configuration to dictionary."""
        return {
            'min_variant_reads': self.min_variant_reads,
            'nucleosome_range': list(self.nucleosome_range) if self.nucleosome_range else None,
            'require_single_variant': self.require_single_variant,
            'min_read_length': self.min_read_length,
            'size_bins': self.size_bins,
            'position_range': list(self.position_range) if self.position_range else None,
            'fdr_threshold': self.fdr_threshold,
            'min_effect_size': self.min_effect_size,
            'output_format': self.output_format,
            'per_variant_output': self.per_variant_output,
            'generate_static_figures': self.generate_static_figures,
            'heatmap_cmap': self.heatmap_cmap,
            'heatmap_vmin': self.heatmap_vmin,
            'heatmap_vmax': self.heatmap_vmax,
            'significance_cmap': self.significance_cmap,
            'significance_vmax': self.significance_vmax,
        }
    
    @classmethod
    def from_dict(cls, config_dict: Dict) -> 'Config':
        """Create Config from dictionary."""
        # Handle tuple conversion
        if 'nucleosome_range' in config_dict and config_dict['nucleosome_range']:
            config_dict['nucleosome_range'] = tuple(config_dict['nucleosome_range'])
        if 'position_range' in config_dict and config_dict['position_range']:
            config_dict['position_range'] = tuple(config_dict['position_range'])
        
        # Filter to only valid fields
        valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
        filtered = {k: v for k, v in config_dict.items() if k in valid_fields}
        
        return cls(**filtered)
    
    def save(self, path: str):
        """Save configuration to YAML file."""
        with open(path, 'w') as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False)
        logger.info(f"Configuration saved to {path}")
    
    def __str__(self) -> str:
        """Pretty-print configuration."""
        lines = ["Configuration:"]
        for key, value in self.to_dict().items():
            lines.append(f"  {key}: {value}")
        return "\n".join(lines)


def load_config(path: str) -> Config:
    """
    Load configuration from a YAML file.
    
    Parameters:
    -----------
    path : str
        Path to YAML configuration file
        
    Returns:
    --------
    Config
        Loaded configuration object
        
    Raises:
    -------
    FileNotFoundError
        If config file doesn't exist
    ValueError
        If config file is invalid
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {path}")
    
    logger.info(f"Loading configuration from {path}")
    
    with open(path) as f:
        config_dict = yaml.safe_load(f)
    
    if config_dict is None:
        config_dict = {}
    
    # Merge with defaults
    merged = DEFAULT_CONFIG.copy()
    merged.update(config_dict)
    
    # Create config object
    config = Config.from_dict(merged)
    
    # Validate
    errors = config.validate()
    if errors:
        raise ValueError(f"Invalid configuration:\n" + "\n".join(f"  - {e}" for e in errors))
    
    return config


def create_default_config_file(path: str):
    """
    Create a default configuration file.
    
    Parameters:
    -----------
    path : str
        Path to write the configuration file
    """
    config = Config()
    config.save(path)
    logger.info(f"Default configuration file created at {path}")


def merge_cli_args(config: Config, args: Any) -> Config:
    """
    Merge command-line arguments into configuration.
    
    CLI arguments override config file values.
    
    Parameters:
    -----------
    config : Config
        Base configuration
    args : argparse.Namespace or dict
        Command-line arguments
        
    Returns:
    --------
    Config
        Updated configuration
    """
    if hasattr(args, '__dict__'):
        args = vars(args)
    
    config_dict = config.to_dict()
    
    # Map CLI argument names to config names (if different)
    cli_to_config = {
        'min_reads': 'min_variant_reads',
        'nuc_min': None,  # Special handling below
        'nuc_max': None,
        'pos_start': None,
        'pos_end': None,
    }
    
    # Handle nucleosome range
    if 'nuc_min' in args and 'nuc_max' in args:
        if args['nuc_min'] is not None and args['nuc_max'] is not None:
            config_dict['nucleosome_range'] = (args['nuc_min'], args['nuc_max'])
    
    # Handle position range
    if 'pos_start' in args and 'pos_end' in args:
        if args['pos_start'] is not None and args['pos_end'] is not None:
            config_dict['position_range'] = (args['pos_start'], args['pos_end'])
    
    # Handle other arguments
    for cli_name, value in args.items():
        if value is None:
            continue
        if cli_name in cli_to_config:
            config_name = cli_to_config[cli_name]
            if config_name:
                config_dict[config_name] = value
        elif cli_name in config_dict:
            config_dict[cli_name] = value
    
    return Config.from_dict(config_dict)


if __name__ == "__main__":
    # Test configuration
    import sys
    
    logging.basicConfig(level=logging.INFO)
    
    # Create default config
    print("Default configuration:")
    config = Config()
    print(config)
    
    # Validate
    errors = config.validate()
    if errors:
        print(f"\nValidation errors: {errors}")
    else:
        print("\nConfiguration is valid!")
    
    # Test saving/loading
    if len(sys.argv) > 1:
        test_path = sys.argv[1]
        config.save(test_path)
        loaded = load_config(test_path)
        print(f"\nLoaded configuration from {test_path}:")
        print(loaded)
