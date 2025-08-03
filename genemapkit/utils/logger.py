"""
Logging Utilities for GeneMapKit

Provides consistent logging configuration across the application.
"""

import logging
import sys
from typing import Optional
import colorlog


def setup_logging(level: str = 'INFO', log_file: Optional[str] = None) -> logging.Logger:
    """
    Set up logging configuration for GeneMapKit.
    
    Args:
        level: Logging level ('DEBUG', 'INFO', 'WARNING', 'ERROR')
        log_file: Optional file path to write logs to
        
    Returns:
        Configured logger instance
    """
    # Create logger
    logger = logging.getLogger('genemapkit')
    logger.setLevel(getattr(logging, level.upper()))
    
    # Clear any existing handlers
    logger.handlers.clear()
    
    # Create formatter
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    # Console handler with colors
    console_handler = colorlog.StreamHandler(sys.stdout)
    console_formatter = colorlog.ColoredFormatter(
        '%(log_color)s%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        log_colors={
            'DEBUG': 'cyan',
            'INFO': 'green',
            'WARNING': 'yellow',
            'ERROR': 'red',
            'CRITICAL': 'red,bg_white',
        }
    )
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    # File handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_formatter = logging.Formatter(log_format, datefmt='%Y-%m-%d %H:%M:%S')
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
    
    return logger


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger instance for a specific module.
    
    Args:
        name: Logger name (usually __name__)
        
    Returns:
        Logger instance
    """
    return logging.getLogger(f'genemapkit.{name}')


class LoggingMixin:
    """Mixin class to add logging capabilities to other classes."""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger = get_logger(self.__class__.__name__)


def log_function_call(func):
    """Decorator to log function calls."""
    def wrapper(*args, **kwargs):
        logger = get_logger(func.__module__)
        logger.debug(f"Calling {func.__name__} with args={args}, kwargs={kwargs}")
        result = func(*args, **kwargs)
        logger.debug(f"{func.__name__} returned: {result}")
        return result
    return wrapper
