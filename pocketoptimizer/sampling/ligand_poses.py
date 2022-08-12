from math import sin, cos
import multiprocessing as mp
from functools import partial
from typing import List, Tuple, Dict, NoReturn
import logging
import os
from copy import deepcopy
import numpy as np
from tqdm.auto import tqdm
from ffevaluation.ffevaluate import FFEvaluate
from moleculekit.molecule import Molecule
from moleculekit.util import uniformRandomRotation
from moleculekit.projections.metricdistance import MetricDistance

logger = logging.getLogger(__name__)


class PoseSampler:
    """
    Class providing methods to create an ensemble of ligand poses inside the binding pocket and to prune generated poses based on a vdw threshold.
    """
    def __init__(self, work_dir: str, ligand: str, forcefield: str = 'amber_ff14SB', sampling_method: str = 'grid',
                 max_poses: int = 10000, filter: str = 'diverse', peptide: bool = False):
        """
        Constructor method.

        Parameters
        ----------
        work_dir: str
            working directory
        ligand: str
            Path to ligand conformer file
        forcefield: str
            Force field that will be used to compute the vdW energy
        sampling_method: str
            Method that will be used to generate poses ('grid', 'random') [default: 'grid']
        max_poses: maximum number of poses to generate [default: 10000]
        filter: filtering method ('diverse', 'distance') [default: 'diverse']
                option: 'diverse' filters diverse poses using the min/max-diversity picker from RDKit
                after pruning of poses and option: 'distance' filters poses based on the distance to scaffold
                before pruning of poses
        peptide: bool
            Whether ligand is a peptide [default: False]
        """
        self.work_dir = work_dir
        self.ligand = ligand
        self.forcefield = forcefield
        self.sampling_method = sampling_method
        self.max_poses = max_poses
        self.filter = filter
        self.peptide = peptide

    @ staticmethod
    def rotate_x(angle: float) -> np.ndarray:
        """
        Calculates Matrix needed to transform ligand coords for rotation of angle degrees around x-axis

        Parameters
        ----------
        angle: angle which will be rotated around the x-axis

        Returns
        -------
        numpy array needed for transforming ligand coordinates
        """
        return np.array([
            [1, 0, 0],
            [0, cos(angle), -sin(angle)],
            [0, sin(angle), cos(angle)]
        ])

    @ staticmethod
    def rotate_y(angle: float) -> np.ndarray:
        """
        Calculates Matrix needed to transform ligand coords for rotation of angle degrees around y-axis

        Parameters
        ----------
        angle: angle which will be rotated around the y-axis

        Returns
        -------
        numpy array needed for transforming ligand coordinates
        """
        return np.array([
            [cos(angle), 0, sin(angle)],
            [0, 1, 0],
            [-sin(angle), 0, cos(angle)]
        ])

    @ staticmethod
    def rotate_z(angle: float) -> np.ndarray:
        """
        Calculates Matrix needed to transform ligand coords for rotation of angle degrees around z-axis

        Parameters
        ----------
        angle: angle which will be rotated around the z-axis

        Returns
        -------
        numpy array needed for transforming ligand coordinates
        """
        return np.array([
            [cos(angle), -sin(angle), 0],
            [sin(angle), cos(angle), 0],
            [0, 0, 1]
        ])

    def rotational_sampling(self, angles: float, steps: float) -> np.ndarray:
        """
        Calculates transformations to rotate ligand in 3 dimensions in a range from -angles to +angles with increment step

        Parameters
        ----------
        angles: maximum and minimum rotation
        steps: increment

        Returns
        -------
        numpy array of matrices needed for the rotational transformations in all 3 dimensions
        """
        sample_angles = np.linspace(-angles, angles, int((angles * 2) / steps + 1))*np.pi / 180
        x_rot = []
        y_rot = []
        z_rot = []
        for angle in sample_angles:
            x_rot.append(self.rotate_x(angle))
            y_rot.append(self.rotate_y(angle))
            z_rot.append(self.rotate_z(angle))
        return np.array(x_rot), np.array(y_rot), np.array(z_rot)

    @ staticmethod
    def rotate(coords: np.ndarray, mat: np.ndarray, center: np.ndarray) -> np.ndarray:
        """
        Calculates new coordinates after rotation

        Parameters
        ----------
        coords: Starting coordinates
        mat: matrix needed for rotation
        center: Center of mass of the molecule

        Returns
        -------
        New coordinates after rotation
        """
        new_coords = coords - center
        return np.dot(new_coords, np.transpose(mat)) + center

    @ staticmethod
    def translate(coords: np.ndarray, vector: np.ndarray) -> np.ndarray:
        """
        Calculates new coordinates after translation
        Parameters
        ----------
        coords: Starting coordinates
        vector: Translation vector

        Returns
        -------
        New coordinates after translation
        """
        return coords + vector

    @ staticmethod
    def translational_sampling(distance: float, steps: float) -> np.ndarray:
        """
        Calculates transfortmations to translate ligand in all 3 dimensions between -distance and +distance with icnrement steps

        Parameters
        ----------
        distance: maximum and minimum distance to translate
        steps: increment

        Returns
        -------
        numpy array of translates coordinates
        """
        sample_range = np.linspace(-distance, distance, int((distance*2)/steps+1))

        x_trans = []
        y_trans = []
        z_trans = []

        for step in sample_range:
            x_trans.append(np.array([step, 0, 0]))
            y_trans.append(np.array([0, step, 0]))
            z_trans.append(np.array([0, 0, step]))

        return np.array(x_trans), np.array(y_trans), np.array(z_trans)

    @ staticmethod
    def rmsd_calc(coords1: np.ndarray, coords2: np.ndarray) -> float:
        """
        Calculates root mean square deviation of two coordinate arrays
        Parameters
        ----------
        coords1: first coordinate array
        coords2: second coordinate array

        Returns
        -------
        RMSD of the two arrays
        """
        dist = coords1 - coords2
        return np.sqrt(np.mean(np.sum(dist * dist, axis=1), axis=0))

    def sample_poses(self, coordinates: np.ndarray, center: np.ndarray, trans: float = 2.0, rot: float = 90, trans_step: float = 0.25, rot_step: float = 15) -> np.ndarray:
        """Calculates all matrix transformations based on the defined rotation and translation.

        First rotates the coordinates along every axis and then translates every created rotation
        along every axis.

        Parameters
        ----------
        coordinates: array
            3D-array containing molecule coordinates.
        center: array
            The center off mass of the molecule.
        trans: float
            Maximum translation distance in Angstrom [default 2.0 Angstrom].
        rot: float
            Maximum rotation of the molecule in degree [default 90 Degree].
        trans_step: float
            Stepsize for the translational sampling [default 0.25 Angstrom].
        rot_step: float
            Stepsize for the rotational sampling [default 15 Degree].

        Returns
        -------
        3D-array containing all the poses.
        """
        rotated_coords = []
        sampled_poses = []

        if rot:
            rot_samples = self.rotational_sampling(rot, rot_step)
            for axis in rot_samples:
                for matrix in axis:
                    rotated_coords.append(self.rotate(coordinates, matrix, center))
        else:
            rotated_coords.append(coordinates)
        if trans:
            trans_samples = self.translational_sampling(trans, trans_step)
            for pose in rotated_coords:
                for axis in trans_samples:
                    for vector in axis:
                        sampled_poses.append(self.translate(pose, vector))
        if not sampled_poses:
            sampled_poses.append(coordinates)

        return np.array(sampled_poses).transpose(1, 2, 0)

    @ staticmethod
    def filter_clashes(pose_id: int, pose_coords: np.ndarray, struc: Molecule, ffev: FFEvaluate) -> np.float:
        """Calculation and filtering procedure.

        Calculates the vdw energy of a ligand pose

        Parameters
        ----------
        pose_id: int
            Pose id
        pose_coords: np.ndarray
            coordinates of all poses
        struc: Molecule class
            moleculekit molecule object.
        ffev: class
            FFevaluate object.

        Returns
        -------
        Vdw energy of a ligand pose
        """
        struc.set('coords', pose_coords[:, :, pose_id], 'segid L')
        energies = ffev.calculateEnergies(struc.coords)
        return energies['vdw']

    @ staticmethod
    def filter_redundant(coordinates: np.ndarray) -> np.ndarray:
        """
        Filters redundant (identical) coordinates that could be created during the sampling.

        I.e. translating from -1 to 1 angstrom in 0.5 angstrom steps produces three poses at
        -1, 0, 1, but 0 in this case would be the starting pose.

        Parameters
        ----------
        coordinates: array
            3-D numpy array with the shape (natoms, xyz, poses)
        Returns
        -------
        coordinates: array
            3-D numoy array that only contains unique poses.
        """
        return np.unique(coordinates, axis=2)

    @ staticmethod
    def filter_distances(struc: Molecule, poses: np.ndarray, min_dist: float = 0.5, max_dist: float = 7.0) -> np.ndarray:
        """
        Filters poses based on distances to the scaffold.

        Parameters
        ----------
        struc: moleculekit Molecule Object
            Molecule Object containing scaffold and a ligand (segid L).
        poses: np.array
            Ligand poses.
        min_dist: float
            Closest distance to the scaffold.
        max_dist: float
            Maximum distance to the scaffold.

        Returns
        -------
        A np.array containing valid poses.
        """
        metric_dist = MetricDistance('protein', 'segid L', pbc=False, groupsel1='all', groupsel2='residue')
        ids = []
        for i in tqdm(range(poses.shape[-1])):
            struc.set('coords', poses[:, :, i], 'segid L')
            if min_dist < metric_dist.project(struc).flatten()[0] < max_dist:
                ids.append(i)

        return poses[..., np.array(ids)]

    def filter_diverse(self, poses: np.ndarray, nsamples: int, diversity: str = 'max', first_picks: List = None) -> np.ndarray:
        """
        Fitlers poses using the minmax diversity picker algorithm in RDkit.

        Parameters
        ----------
        poses: np.array
            Ligand poses.
        nsamples: int
            Number of poses to sample.
        diversity: str
            Maximum or minimum diversity. [default: 'max']
        first_picks: list
            List of Ids that are the starting point for the diversity picker.

        Returns
        -------
        A np.array containing valid poses.
        """
        from rdkit.SimDivFilters import rdSimDivPickers
        from itertools import combinations

        nposes = poses.shape[-1]

        lower_matrix = list(combinations(np.arange(nposes), 2))
        rmsds = np.zeros((nposes, nposes))
        for i, j in tqdm(lower_matrix):
            rmsds[j, i] = self.rmsd_calc(poses[..., i], poses[..., j])

        if diversity == 'max':
            rmsds = rmsds[np.tril_indices(nposes, k=1)]
        elif diversity == 'min':
            reverted_rmsd = -rmsds - np.min(rmsds)
            rmsds = reverted_rmsd[np.tril_indices(nposes, k=1)]
        else:
            logger.error(f'{diversity} method not defined. Try "max" or "min".')
            raise ValueError(f'{diversity} method not defined. Try "max" or "min".')

        minmax_picker = rdSimDivPickers.MaxMinPicker()
        fp = list(first_picks) if first_picks is not None else []
        idx = list(minmax_picker.Pick(rmsds, nposes, nsamples, firstPicks=fp))
        logger.info(f'Filtered to {len(idx)} poses using minmax diversity picker.')
        return np.array(idx)

    def grid_sampling(self, mol: Molecule, grid: Dict[str, List[float]]) -> np.ndarray:
        """
        Samples all ligand conformers along the defined grid.

        Parameters
        ----------
        mol: moleculekit Molecule Object
            The ligand.
        grid: dict
            Dictionary containing grid parameters for sampling poses.
            Has the following structure:
             {'trans':[DISTANCE, STEPSIZE],
              'rot':[ANGLE, STEPSIZE]}

        Returns
        -------
        A np.array of all created poses.
        """
        poses = []
        logger.info('Generate possible poses within the defined grid.')
        for conf in tqdm(range(mol.coords.shape[-1]), desc='Ligand Conformers'):
            coords = mol.coords[:, :, conf]
            mol_center = np.mean(coords, axis=0)
            poses.append(self.sample_poses(
                coordinates=coords,
                center=mol_center,
                trans=grid['trans'][0],
                trans_step=grid['trans'][1],
                rot=grid['rot'][0],
                rot_step=grid['rot'][1]
            ))

        return self.filter_redundant(np.dstack(poses))

    def random_sampling(self, mol: Molecule, nposes: int) -> np.ndarray:
        """
        Random sampling procedure for the ligand.

        Parameters
        ----------
        mol: moleculekit Molecule Object
        nposes: int
            Number of poses to create.

        Returns
        -------
        A np.array of all created poses.
        """
        poses = []
        while len(poses) < nposes:
            mol2 = deepcopy(mol)
            for conf in range(mol2.coords.shape[-1]):
                mol2.dropFrames(keep=0)
                mol2.set('coords', mol.coords[:, :, conf])
                mol_center = np.mean(mol2.get('coords'), axis=0)
                rotation_matrix = uniformRandomRotation()
                mol2.rotateBy(rotation_matrix, center=mol_center)
                poses.append(mol2.coords)
        return self.filter_redundant(np.dstack(poses))

    def create_ligand_ensemble(self, grid: Dict[str, List[float]] = None, vdw_filter_thresh: float = 100.0, ncpus: int = 1) -> NoReturn:
        """Creates a ligand pose ensemble by transforming the ligand along every axis.

        Every ligand conformer provided will be translated and rotated along x, y, z
        Afterwards the created poses will be filtered according to a vdW threshold value
        to secure that no pose is clashing with the backbone. The pose sampling procedure
        is performed in a glycin or alanine scaffold.
        After the sampling is finished the created poses will be written into a .xtc
        trajectory file.

        Parameters
        ----------
        grid: dict
            Contains translational and rotational steps.
        vdw_filter_thresh: float
            Threshold value for filtering poses [default: 100 kcal/mol]
        ncpus: int
            Number of CPUs to be used. [default: 1]
        """
        from pocketoptimizer.utility.utils import load_ff_parameters, write_energies, calculate_chunks

        os.makedirs(os.path.join(self.work_dir, 'ligand', self.forcefield, 'poses'), exist_ok=True)

        logger.info('Start ligand pose sampling procedure.')
        struc, prm = load_ff_parameters(structure_path=os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_params', 'ligand_sampling_pocket'), forcefield=self.forcefield)

        try:
            ligand_confs = Molecule(self.ligand)
            # Dont append native conf for peptides as it is already included
            if not self.peptide:
                native_conf = np.expand_dims(struc.get('coords', sel='segid L'), axis=2)
                ligand_confs.coords = np.dstack((ligand_confs.coords, native_conf))
        except (IndexError, FileNotFoundError):
            # If ligand confs doesn't exist or doesn't contain conformers
            logger.error(f'{self.ligand} does not exist or does not contain conformers.')
            raise FileNotFoundError(f'{self.ligand} does not exist or does not contain conformers.')

        # Atom order needs to be the same in both files
        if (struc.get('name', 'segid L') != ligand_confs.get('name')).all():
            logger.error('Unequal atom order between generated ligand conformers and protonated ligand.')
            raise RuntimeError('Unequal atom order between generated ligand conformers and protonated ligand.')

        if self.sampling_method == 'grid':
            poses = self.grid_sampling(ligand_confs, grid)
        elif self.sampling_method == 'random':
            poses = self.random_sampling(ligand_confs, nposes=self.max_poses)
        else:
            logger.error(f'{self.sampling_method} is not a valid method. Use method grid or random.')
            raise ValueError(f'{self.sampling_method} is not a valid method. Use method grid or random.')

        if self.filter == 'distance':
            poses = self.filter_distances(struc, poses)

        nposes = poses.shape[-1]
        logger.info(f'Created possible {nposes} poses.')

        logger.info('Start filtering poses.')
        logger.info(f'Using {ncpus} CPUs for multiprocessing.')

        ffev = FFEvaluate(struc, prm, betweensets=('protein', 'segid L'))
        chunksize = calculate_chunks(nposes=nposes, ncpus=ncpus)

        energies = np.ndarray(nposes)
        with tqdm(total=nposes, desc='Filter Poses') as pbar:
            with mp.Pool(processes=ncpus) as pool:
                    for pose_id, energy in enumerate(pool.imap(
                            partial(self.filter_clashes,
                                    pose_coords=poses,
                                    struc=struc.copy(),
                                    ffev=ffev
                                    ), np.arange(nposes), chunksize=chunksize)):
                        energies[pose_id] = energy
                        pbar.update()

        val_ids = [val_id[0] for val_id in np.argwhere(energies <= min(energies) + vdw_filter_thresh)]
        logger.info(f'Calculated poses within threshold: {len(val_ids)}.')

        ligand_pose_energy_file = os.path.join(self.work_dir, 'ligand', self.forcefield, 'poses', 'ligand_poses.csv')
        write_energies(outpath=ligand_pose_energy_file,
                       energies=energies,
                       energy_terms=['VdW'],
                       name_a='ligand_pose',
                       nconfs_a=nposes)

        val_ids = np.hstack(val_ids)
        if val_ids.shape[0] > self.max_poses:
            if self.filter == 'diverse':
                val_ids = self.filter_diverse(poses[:, :, val_ids], nsamples=self.max_poses)
        ligand_confs.coords = poses[:, :, val_ids]

        # Save as trajectory
        ligand_confs.write(os.path.join(self.work_dir, 'ligand', self.forcefield, 'poses', f'ligand_poses.pdb'))
        ligand_confs.write(os.path.join(self.work_dir, 'ligand', self.forcefield, 'poses', f'ligand_poses.xtc'))
        logger.info('Pose sampling procedure was successful.')
