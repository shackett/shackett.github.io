import logging
import sys

def setup_logging():
    """Set up logging configuration with debug level"""
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(levelname)s: %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
