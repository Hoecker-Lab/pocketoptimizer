import os
import tempfile as tf
from shutil import copyfile, rmtree
import subprocess
from tqdm.auto import tqdm
from typing import Dict, Union, NoReturn
import logging

from pocketoptimizer.utility.index_mapper import IndexMapper

logger = logging.getLogger(__name__)


def run_sontag_solver(sol: int, sontag_params: Dict[str, Union[str, int, float]], output_dir: str, \
                      ligand_index: int, pair_count: int, penalty_energy: float) -> NoReturn:
    """
    Calculates one design solution

    Parameters
    ----------
    sol: int
        Index of the solution to be calculated
    sontag_params: dict
        Keys are names of sontag command line parameters,
        values their values (key "bin" points to path of sontag solver binary)
    output_dir: str
        Path where solution files are written to
    ligand_index: int
        Index of design position "ligand" (among the alphabetically sorted design positions)
    pair_count: int
        Number of design position pairs
    penalty_energy: float
        Dummy energy to be used to exclude a ligand pose from being part of the solution
    """
    # Syntax: algo_triplet.exe niter niter_later nclus_to_add obj_del_thr
    #                            int_gap_thr graph_type (num_rows num_cols)
    command = [
        sontag_params['bin'],
        str(sontag_params['niter']),
        str(sontag_params['niter_later']),
        str(sontag_params['nclust_to_add']),
        str(sontag_params['obj_del_thr']),
        str(sontag_params['int_gap_thr']),
        '0'
    ]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    # check if result was produced and store it
    result_filename = 'res.txt'

    if not os.path.isfile(result_filename):
        logger.error(f'Sontag solver failed with the following exception: {stderr.decode("ascii")}.')
        raise RuntimeError(f'Sontag solver failed with the following exception: {stderr.decode("ascii")}.')

    copyfile(result_filename, os.path.join(output_dir, 'res{sol:0=2d}.txt'))
    with open(result_filename) as resfile:
        solution_line = resfile.readlines()[0]

    with open(os.path.join(output_dir, 'all_solutions.txt'), 'a') as out:
        out.write(solution_line)

    pose_id = int(solution_line.split()[ligand_index])
    os.rename(result_filename, f'res{sol:0=2d}.txt')
    os.rename('lambdas.txt', f'lambdas{sol:0=2d}.txt')
    with open('lambdas.txt', 'w') as lambda_out:
        for i, line in enumerate(open(f'lambdas{sol:0=2d}.txt')):
            # would be the ligand self energy line
            if i == pair_count + ligand_index:
                sp = line.split()
                line = ' '.join(sp[:pose_id]) \
                       + ' %g ' % (-1 * penalty_energy) \
                       + ' '.join(sp[pose_id + 1:]) + '\n'
            lambda_out.write(line)
            lambda_out.flush()


def exclude_solutions(exclude: str, ligand_index: int, pair_count: int, penalty_energy: float) -> NoReturn:
    """
    Exclude solutions

    Parameters
    ----------
    exclude: str
        Path to file containing pose Ids to exclude
    ligand_index: int
        Index of design position "ligand" (among the alphabetically sorted design positions)
    pair_count: int
        Number of design position pairs
    penalty_energy: float
        Dummy energy to be used to exclude a ligand pose from being part of the solution
    """
    os.rename('lambdas.txt', 'lambdas_pre_exclusion.txt')
    with open('lambdas.txt', 'w') as lambda_out:
        poses_to_exclude = []
        with open(exclude) as file_handle:
            for line in file_handle:
                if line.isspace() or not line:
                    continue
                pose_id = int(line.split()[ligand_index])
                poses_to_exclude.append(pose_id)
            with open('lambdas_pre_exclusion.txt') as file_handle_pre:
                for i, line in enumerate(file_handle_pre):
                    # whould be the ligand self energy line
                    if i == pair_count + ligand_index:
                        sp = line.split()
                        for pi in poses_to_exclude:
                            sp[pi] = '%g' % (-1 * penalty_energy)
                        line = ' '.join(sp) + '\n'
                    lambda_out.write(line)
                    lambda_out.flush()


def calculate_design_solutions(solver_bin: str, temp_dir: str, out_path: str, num_solutions: int = 10,
                               penalty_energy: float = 1e10, exclude: str = None, keep_tmp: bool = False) -> NoReturn:
    """
    Iteratively submit design jobs with the Sontag solver, allowing each ligand pose to appear only once in a solution.

    Parameters
    ----------
    solver_path: str
        Path to solver binary
    temp_dir: str
        Path to temporary directory storing the sontag solutions
    out_path: str
        Output directory to write the solutions
    num_solutions: int
        Number of solutions to calculate
    penalty_energy: float
        The dummy energy to be used to exclude a ligand pose from being part of the solution
    exclude: str
        Path to file containing pose Ids to exclude
    keep_tmp: bool
        Whether the temporary directory should be removed or kept in the end [default: False]

    """
    logger.info('Calculating Solutions.')
    index_mapper = IndexMapper.from_index_file(os.path.join(out_path, 'index.dat'))

    if index_mapper.get_conf_count_for_pos('ligand') < num_solutions:
        num_solutions = index_mapper.get_conf_count_for_pos('ligand')
        logger.warning(f'Requested number of solutions higher than number of possible ligand poses. '
              f'Only {num_solutions} can be calculated')

    ligand_index = index_mapper.get_pos_index('ligand')

    # calculate the number of position/pose pairs in the design problem
    pair_count = index_mapper.get_pos_count()
    pair_count = int(((pair_count * pair_count) - pair_count) / 2)

    tmp_dir = tf.mkdtemp(prefix='calcSolutionsWithSontag_', dir=temp_dir)
    os.chdir(tmp_dir)

    # since the solutions.txt file is opened in append mode, if it exists, it should be deleted,
    # otherwise different results might accumulate if called multiple times
    try:
        os.remove(os.path.join(out_path, 'solutions', 'all_solutions.txt'))
    except FileNotFoundError:
        pass

    input_files = ['regions.txt', 'intersects.txt', 'var_sizes.txt',
                   'region_intersects.txt', 'lambdas.txt']
    for in_file in input_files:
        assert os.path.exists(os.path.join(out_path, in_file)), f'Missing input file: {in_file}'
        copyfile(os.path.join(out_path, in_file), in_file)

    if exclude:
        exclude_solutions(exclude, ligand_index, pair_count, penalty_energy)

    sontag_params = {'bin': solver_bin,
                     'niter': 1000,
                     'niter_later': 20,
                     'nclust_to_add': 5,
                     'obj_del_thr': 2e-4,
                     'int_gap_thr': 2e-5}

    for sol in tqdm(range(num_solutions), desc='Solutions'):
        run_sontag_solver(
            sol=sol,
            sontag_params=sontag_params,
            output_dir=os.path.join(out_path, 'solutions'),
            ligand_index=ligand_index,
            pair_count=pair_count,
            penalty_energy=penalty_energy
        )

    if not keep_tmp:
        if os.path.isdir(tmp_dir):
            rmtree(tmp_dir)