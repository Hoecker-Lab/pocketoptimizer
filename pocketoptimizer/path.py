import os
import logging

logger = logging.getLogger(__name__)


def path(bin_file: str = None, rotamers: bool = False) -> str:
    """
    Returns the path name to the PocketOptimizer root directory or
    a specified bin file.


    Parameters
    ----------
    bin_file: str
        Name of a file contained in pocketoptimizer/bin.
    rotamers: bool
        Whether to return the rotamer directory
    Returns
    -------
    path: str
        Path to a directory or file.

    """
    po_path = os.path.split(os.path.realpath(__file__))[0]
    if bin_file:
        bin_path = os.path.join(po_path, 'bin', bin_file)
        if os.path.exists(bin_path):
            return bin_path
        else:
            logger.error('File not found.')
            raise FileNotFoundError('File not found.')
    elif rotamers:
        return os.path.join(os.path.split(po_path)[0], 'rotamers')
    else:
        return po_path
