"""Utility functions for spatial biopsy analysis."""

import os
import logging
import psutil


def log_memory_usage():
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 1024 / 1024
    logging.info(f"Memory usage: {mem:.2f} MB")
