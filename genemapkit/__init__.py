"""
GeneMapKit Package Initialization
"""

__version__ = "1.0.0"
__author__ = "GeneMapKit Development Team"
__description__ = "Comprehensive Gene ID Mapping Toolkit"

from .core.converter import GeneConverter, MAPPING_POLICIES, MappingResult
from .core.downloader import DatabaseDownloader
from .core.external_validation import compare_mygene
from .core.external_benchmark import benchmark_external_services
from .core.publication_benchmark import run_publication_benchmark
from .core.path_benchmark import benchmark_conversion_paths

__all__ = [
    'GeneConverter',
    'MappingResult',
    'MAPPING_POLICIES',
    'DatabaseDownloader',
    'compare_mygene',
    'benchmark_external_services',
    'run_publication_benchmark',
    'benchmark_conversion_paths',
]
