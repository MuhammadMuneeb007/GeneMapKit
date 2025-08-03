"""
GeneMapKit Package Initialization
"""

__version__ = "1.0.0"
__author__ = "GeneMapKit Development Team"
__description__ = "Comprehensive Gene ID Mapping Toolkit"

from .core.converter import GeneConverter
from .core.downloader import DatabaseDownloader

__all__ = ['GeneConverter', 'DatabaseDownloader']
