import os
import logging

logger = logging.getLogger(__name__)


def path(data_dir: str = None, bin_file: str = None) -> str:
    """
    Returns the path name to the PocketOptimizer root directory or
    a specified data directory or bin file.


    Parameters
    ----------
    data_dir: str
        Name of a directory contained in pocketoptimizer/data.
    bin_file: str
        Name of a file contained in pocketoptimizer/bin.

    Returns
    -------
    path: str
        Path to a directory or file.

    """
    po_path = os.path.split(os.path.realpath(__file__))[0]
    if data_dir:
        dpath = os.path.join(po_path, 'data', data_dir)
        if os.path.exists(dpath):
            return dpath
        else:
            logger.error('Data path not found.')
            raise FileNotFoundError('Data path not found.')
    elif bin_file:
        bin_path = os.path.join(po_path, 'bin', bin_file)
        if os.path.exists(bin_path):
            return bin_path
        else:
            logger.error('File not found.')
            raise FileNotFoundError('File not found.')
    else:
        return po_path
