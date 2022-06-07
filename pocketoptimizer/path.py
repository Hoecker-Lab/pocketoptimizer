import os
import shutil
import glob
import logging
from typing import NoReturn

logging.root.handlers = []
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - [%(levelname)s] - %(message)s",
    handlers=[
        logging.FileHandler(os.environ.get('POCKETOPTIMIZER_LOGFILE')),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger('pocketoptimizer.path')


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


def data_path(data_dir: str, target: str) -> str:
    """
    Returns the path to subdirectories inside the data directory.

    Parameters
    ----------
    data_dir: str
        Valid options are 'benchmark', 'ligands', 'structures', 'tests'.
    target: str
        Subdirectory inside the data_dir.

    Returns
    -------
    dir: str
        Path to the directory

    """

    valid_paths = ['benchmark', 'ligands', 'structures', 'tests']
    if data_dir not in valid_paths:
        logger.error(f'data_dir must be one of {valid_paths}.')
        raise ValueError(f'data_dir must be one of {valid_paths}.')

    bench_path = os.path.join(path(data_dir=data_dir), target)
    if os.path.exists(bench_path):
        return bench_path
    else:
        logger.error('Data path not found.')
        raise FileNotFoundError('Data path not found.')


def copy_data(source: str, destination: str) -> NoReturn:
    """
    Copies a project directory

    Parameters
    ----------
    source: str
        Path to original source directory
    destination: str
        Path to the new copy
    """

    os.makedirs(destination, exist_ok=True)
    os.makedirs(os.path.join(destination, 'scaffold'), exist_ok=True)
    os.makedirs(os.path.join(destination, 'ligand'), exist_ok=True)
    lig_files = glob.glob(os.path.join(source, 'ligand', '*'))
    scaff_files = glob.glob(os.path.join(source, 'scaffold', '*'))
    for i in lig_files:
        fname = os.path.basename(i)
        shutil.copy(i, os.path.join(destination, 'ligand', fname))
    for i in scaff_files:
        fname = os.path.basename(i)
        shutil.copy(i, os.path.join(destination, 'scaffold', fname))
