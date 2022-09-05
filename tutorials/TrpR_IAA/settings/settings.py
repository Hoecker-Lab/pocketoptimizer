OBABEL_BIN = '/agh/projects/jakob/miniconda3/envs/pocketoptimizer/bin/obabel'
MATCH_PERL = '/agh/projects/jakob/pocketoptimizer/pocketoptimizer/../MATCH_RELEASE/MATCH/scripts/MATCH.pl'
ANTECHAMBER_BIN = '/agh/projects/jakob/miniconda3/envs/pocketoptimizer/bin/antechamber'
PARMCHK2_BIN = '/agh/projects/jakob/miniconda3/envs/pocketoptimizer/bin/parmchk2'
TLEAP_BIN = '/agh/projects/jakob/miniconda3/envs/pocketoptimizer/bin/tleap'
PSFGEN_BIN = '/agh/projects/jakob/pocketoptimizer/pocketoptimizer/bin/psfgen'
SMINA_BIN = '/agh/projects/jakob/miniconda3/envs/pocketoptimizer/bin/smina'
SOLVER_BIN = '/agh/projects/jakob/pocketoptimizer/pocketoptimizer/bin/sontag_solver'
POCKETOPTIMIZER_LOGFILE = '/agh/projects/jakob/pocketoptimizer/tutorials/TrpR_IAA/pocketoptimizer.log'
TMP_DIR = '/tmp/user/1702035'

class Settings:

    def __init__(self):
        self.OBABEL_BIN = OBABEL_BIN
        self.MATCH_PERL = MATCH_PERL
        self.ANTECHAMBER_BIN = ANTECHAMBER_BIN
        self.PARMCHK2_BIN = PARMCHK2_BIN
        self.TLEAP_BIN = TLEAP_BIN
        self.PSFGEN_BIN = PSFGEN_BIN
        self.SMINA_BIN = SMINA_BIN
        self.SOLVER_BIN = SOLVER_BIN
        self.POCKETOPTIMIZER_LOGFILE = POCKETOPTIMIZER_LOGFILE
        self.TMP_DIR = TMP_DIR