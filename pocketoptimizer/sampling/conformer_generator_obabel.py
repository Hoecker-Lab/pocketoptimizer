import os
import subprocess
import logging
from typing import NoReturn

logger = logging.getLogger(__name__)


def merge_confs(infile: str, outfile: str) -> NoReturn:
    """
    Merges all the output conformers of OBabel into the standard
    .pdb Model format.

    Parameters
    ----------
    infile: str
        Input file.
    outfile: str
        Output file.
    """
    with open(infile, 'r') as tmpfile, open(outfile, 'w') as output:
        counter = 1
        for line in tmpfile.readlines():
            if line.startswith('COMPND'):
                output.write(f'MODEL        {counter}\n')
                counter += 1
            elif line.startswith('END'):
                output.write(f'ENDMDL\n')
            else:
                output.write(line)
        output.write('END\n')
    os.remove(infile)
    assert os.path.isfile(outfile)
    logger.info(f'Generated {counter-1} conformers.')


def generate_genetic(obabel_path: str, infile: str, outfile: str, nconfs: int = 100, score: str = 'rmsd') -> NoReturn:
    """
    Uses a genetic algorithm in OBabel, that finds a globally optimum solution to a multiparameter
    problem from a population of conformers after a series of generations

    If you want all the OBabel output logs you can delete the stdout and stderr parts. Can be helpful for
    debugging.

    Parameters
    ----------
    obabel_path: str
            path to obabel binary
    infile: str
        Input file.
    outfile: str
        Output file.
    nconfs: number of conformers to sample [default: 100]
    score: arriving at an optimal solution in terms of either rmsd diversity or energy [default: 'rmsd']

    """
    f_ext = os.path.splitext(infile)[-1][1:]
    conf_command = [obabel_path, f'-i{f_ext}', infile, '-opdb', '-O', outfile, '--conformer', '--nconf', str(nconfs),
                    '--score', score, '--writeconformers']

    process = subprocess.Popen(conf_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        logger.error(f'Obabel failed with the following error: {stderr.decode("ascii")}')
        raise RuntimeError(f'Obabel failed with the following error: {stderr.decode("ascii")}')
    assert os.path.isfile(outfile)


def generate_confab(obabel_path: str, infile: str, outfile: str, nconfs: int = 100, rcutoff: float = 0.5, ecutoff: float = 50.0) -> NoReturn:
    """
    Uses the confab procedure in OBabel, that systematically generates all diverse low-energy conformers

    If you want all the OBabel output logs you can delete the stdout and stderr parts. Can be helpful for
    debugging.

    Parameters
    ----------
        obabel_path: str
                path to obabel binary
    infile: str
        Input file.
    outfile: str
        Output file.
    nconfs: number of conformers to sample [default: 100]
    rcutoff: RMSD cutoff [default: 0.5 Angstrom]
    ecutoff: energy cutoff [default: 50 kcal/mol]

    """
    f_ext = os.path.splitext(infile)[-1][1:]
    conf_command = [obabel_path, f'-i{f_ext}', infile, '-opdb', '-O', outfile, '--confab', '--conf', str(nconfs),
                    '--rcutoff', str(rcutoff), '--ecutoff', str(ecutoff)]

    process = subprocess.Popen(conf_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        logger.error(f'Obabel failed with the following error: {stderr.decode("ascii")}')
        raise RuntimeError(f'Obabel failed with the following error: {stderr.decode("ascii")}')
    assert os.path.isfile(outfile)


def conformer_generator(obabel_path: str, infile: str, conf_file_name: str, method: str = 'genetic', nconfs: int = 100, score: str = 'rmsd', rcutoff: float = 0.5,
                        ecutoff: float = 50.0) -> NoReturn:
    """Conformer generation using procedures implemented in Obabel.

    Possible methods are 'genetic' or 'confab'. The output conformers will be superimposed and
    saved into a .pdb file.

    Parameters
    ----------
    obabel_path: str
                path to obabel binary
    infile: str
        Path to the input ligand
    conf_file_name: str
        Output conformer file name
    method: str
        Method that will be used to generate conformers. Options are 'genetic' or 'confab'[default: 'genetic].
    nconfs: int
        Maximum number of conformers to sample. [default: 100]
    score: str
        Scoring option to decide how to filter the generated conformers. Only used with the 'rdkit' method. [default: 'rmsd']
    rcutoff: float
        RMSD cutoff [default: 0.5 Angstrom]. Only used together with the 'confab' method.
    ecutoff: float
        Energy cutoff [default 50.0 kcal/mol]. Only used together with the 'confab' method.

    """

    tmp_outfile = os.path.join(os.path.dirname(conf_file_name), f'tmp_conf.pdb')

    if os.path.isfile(conf_file_name):
        logger.info('Conformers are already sampled.')
        return

    logger.info('Starting ligand conformer generation using obabel.')
    logger.info(f'Selected Method: {method}.')

    if method == 'genetic':
        generate_genetic(obabel_path=obabel_path, infile=infile, outfile=tmp_outfile, nconfs=nconfs, score=score)
    elif method == 'confab':
        generate_confab(obabel_path=obabel_path, infile=infile, outfile=tmp_outfile, nconfs=nconfs, rcutoff=rcutoff, ecutoff=ecutoff)
    else:
        logger.error('Method does not exist.')
        raise ValueError('Method does not exist.')
    merge_confs(infile=tmp_outfile, outfile=conf_file_name)
    logger.info('Conformer sampling was successful.')