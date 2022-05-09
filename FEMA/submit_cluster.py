#!/bin/python
#cython: language_level=3
"""Runs FEMA permutations on a distributed way on the specified cluster
Example use:
    python3 submit_cluster.py -h (get help)
    python3 submit_cluster.py -dm /path/to/design_matrix.txt -dt vertex -mod smri -img area -cols "[11,12]" -o /path/to/outputdir
"""

import argparse
import subprocess
from pathlib import Path
import logging
import re

IMAGING_TYPE_DICT = {
    "thickness": {
        "3.0": "thickness-sm{sm}",
        "4.0": "thickness_ic5_sm{sm}"
    },
    "area": {
        "3.0": "area-sm{sm}",
        "4.0": "area_ic5_sm{sm}"
    }
}

class Analysis:
    """Analysis class"""
    def __init__(self,
        design_matrix:Path,
        datatype:str,
        modality:str,
        fstem_imaging:str,
        colsinterest:str, 
        contrasts:list,
        ranknorm:int,
        random_effects:str,
        tfce:int,
        mediation:int,
        permutation_type:str,
        ico:int,
        smoothing_factor:int,
        data_release:str
    ):
        self.design_matrix = design_matrix.absolute()
        self.datatype = datatype
        self.modality = modality
        self.colsinterest = colsinterest           # Columns in design matrix to loop over to calculate TFCE - selecting columns of interest improves efficiency - is `isempty(tfce_cols)` FEMA_wrapper will NOT run TFCE
        self.contrasts = contrasts                 # Contrasts relate to columns in design matrix e.g. [1 -1] will take the difference between cols 1 and 2 in your design matrix (X)
        self.ranknorm = ranknorm                   # Rank normalizes dependent variables (Y) (default = 0)
        self.random_effects = random_effects       # Random effects to include: family, subject, error
        self.tfce = tfce
        self.mediation = mediation                 # If wanting to use outputs for a mediation analysis set mediation=1 - ensures same resampling scheme used for each model in fname_design
        self.permutation_type = permutation_type   # Default resampling method is null wild-bootstrap - to run mediation analysis need to use non-null wild-bootstrap ('wildboostrap-nn')
        self.ico = ico
        self.smoothing_factor = smoothing_factor
        self.data_release = data_release
        # If vertex & smri analysis, construct fstem_imaging, otherwise, use as is.
        if self.datatype == 'vertex' and self.modality == 'smri':
            self.fstem_imaging = IMAGING_TYPE_DICT[fstem_imaging][self.data_release].format(sm=self.smoothing_factor)
        else:
            self.fstem_imaging = fstem_imaging

    def recap(self):
        """Create string with main parameters"""
        to_keep = [
            f'rel-{self.data_release}',
            f'dt-{self.datatype}',
            f'mod-{self.modality}',
            f'img-{self.fstem_imaging}',
        ]
        return '_'.join(to_keep)
    
    def get_genetics(self, abcd_sync:Path) -> Path:
        """Get path to genetic pihat"""
        genp = abcd_sync / self.data_release / f'genomics/ABCD_rel{self.data_release}_pihat.mat'
        if not genp.exists():
            genp = abcd_sync / '3.0' / 'genomics/ABCD_rel3.0_pihat.mat'
            # raise NotFoundErr(genp)
        return genp


class Job:
    """Job object"""
    def __init__(self, name:str, code_dir:Path, command:str) -> None:
        self.name = name
        self.code_dir = code_dir
        self.command = self.matlab_call(command)


    def matlab_call(self, matlab_command:str):
        """Checkproof a matlab command before bash submission"""
        matlab_command = matlab_command.replace(" ","")
        call_command = f"{self.code_dir}/FEMA/submit_cluster_do {self.code_dir} {matlab_command};"
        return call_command

class Queue:
    """Queue object"""
    RESOURCES_DICT = {
        "SGE": "-l mem_free={ram}G -l h_rt={h_rt}",
        "SLURM": "--mem-per-cpu={ram}G --time={h_rt}"
    }

    def __init__(self, cluster:str,mem_free:float,run_time:str) -> None:
        # cluster_type (str): between SLURM, SGE, LOCAL.
        self._cluster = cluster
        self.resources = self.RESOURCES_DICT[cluster].format(ram=mem_free, h_rt=run_time)

    @property
    def cluster(self):
        """`cluster` property"""
        return self._cluster

    @cluster.setter
    def cluster(self, value):
        accepted_clusters = ['SLURM', 'SGE', 'LOCAL']
        if value not in accepted_clusters:
            raise ValueError(f"Accepted clusters: {accepted_clusters}")
        self._cluster = value

    def check_submission(self, output: subprocess.Popen) -> str:
        """Check output of submission call for job_id

        Args:
            output (subprocess.Popen): output of submission call via subprocess.run

        Raises:
            Exception: the job submission failed

        Returns:
            str: job id of the submitted job
        """

        if self.cluster == 'SLURM':
            job_id = output.stdout.rstrip()
            if job_id!= '':
                return job_id
            else:
                raise Exception(
                    f'Error in job submission. '
                    f'Job name: {output.args[2]}\n'
                    f'{output.stderr}')
        elif self.cluster == 'SGE':
            regex_pattern = re.compile(r'(\d+) \("(.*)"\)')
            try:
                job_id, job_name = re.findall(regex_pattern, output.stdout.rstrip())[0]
            except Exception as exc:
                logging.error(f'Output stderr: \n{output.stderr}')
                logging.error(f'Script error: {exc}')
                raise RuntimeError
            return job_id
        else:
            return 'LOCAL'

    def submit_job(self,
                   job:Job,
                   log_dir: Path,
                   dependencies:list,
                   loop_argument: str='',
                   ) -> str:
        """Submit a given job for computation

        Args:
            dependencies (List[str]): other job ids upon which this one depends
            log_dir (Path): directory to store the job output&error logs
            loop_argument (str, optional): if the job is submitted for several cases. Defaults to ''.

        Raises:
            ValueError: if unknown cluster type

        Returns:
            str: successfully submitted job id
        """

        if self.cluster == 'SLURM':
            if len(dependencies) == 0:
                dependencies_arg = ''
            else:
                dependencies_arg = ' --dependency afterok:' + ':'.join(dependencies)
            submit_command = (
                f"sbatch --parsable"
                f" --job-name={job.name}"
                f"{dependencies_arg}"
                f" {self.resources}"
                f" -o {log_dir}/{job.name}.o%j"
                f" {job.command}"
                f" {loop_argument}"
            )
        elif self.cluster == 'SGE':
            if len(dependencies) == 0:
                dependencies_arg = ''
            else:
                dependencies_arg = ' -hold_jid ' + ','.join(dependencies)
            submit_command = (
                f"qsub -o {log_dir} -j y"
                f" -N {job.name}"
                f"{dependencies_arg}"
                f" {self.resources}"
                f" {job.command}"
                f" {loop_argument}"
            )
        elif self.cluster == 'LOCAL':
            submit_command = f'./{job.command} {loop_argument} '
        else:
            raise ValueError(f'Unknown cluster type: {self.cluster}')

        log_dir.mkdir(exist_ok=True, parents=True)
        logging.debug(f"submit_command={submit_command}")
        output = subprocess.run(submit_command.split(' '), capture_output=True, text=True)
        if self.cluster == 'LOCAL':
            log_path = log_dir / f'{job.name}_local.00.log'
            index = 0
            while log_path.exists():
                index += 1
                log_path = log_path.with_name(
                    log_path.stem.split(f'.{index-1:02}')[0] + f'.{index:02}.log')
            with open(log_path, 'w+') as log_f:
                log_f.write(output.stderr)
                log_f.write(output.stdout)

        job_id = self.check_submission(output)
        logging.info(f"Correctly submitted job_id={job_id}, jobname={job.name}")
        return job_id

class Matlab:
    """Wrapper for matlab command calls"""
    def __init__(self, matlab_function, *args, **kwargs) -> None:
        self.function = matlab_function
        self.args = args
        self.kwargs = kwargs
        self.command = self.create_command()

    def create_command(self):
        """Create matlab command as string"""
        def add_value(mystr:str, value) -> str:
            if isinstance(value, (int, float)) or str(value).startswith(('{', '[')):
                mystr += str(f"{value},")
            else:
                mystr += str(f"'{value}',")
            return mystr
        mystr = f"{self.function}("
        # Add main arguments
        for value in self.args:
            mystr = add_value(mystr, value)
        # Add optional arguments
        for key, value in self.kwargs.items():
            mystr += str(f"'{key}',")
            mystr = add_value(mystr, value)
        return mystr.rstrip(',') + ');'

def check_repos(code_dir:Path):
    """Checking if the repos needed are present"""
    repos = ['FEMA', 'cmig_tools_utils'] #'PALM'
    if all((code_dir / repo).is_dir() for repo in repos):
        return
    else:
        raise FileNotFoundError(f'Missing repos in {code_dir} from {*repos,}.')

def chunks(lst, n):
    """Yield successive n-sized chunks from lst"""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def matlab_list(lst:list) -> str:
    """Return Matlab list as string from Python list"""
    lst = list(map(str,lst))
    return str(lst).replace("[", "{").replace("]", "}").replace(" ","")

def set_imaging_dir(abcd_sync, analysis_args:Analysis) -> Path:
    """Setting the correct path for the imaging_dir"""
    if analysis_args.datatype == 'voxel':
        if analysis_args.data_release == '3.0':
            datatype = 'voxelwise/ABCD1_cor10/volmats_subsamp'
        elif analysis_args.data_release == '4.0':
            datatype = 'voxel/ABCD2_cor10'
    else:
        datatype = analysis_args.datatype
    return abcd_sync / analysis_args.data_release / 'imaging_concat' / datatype / analysis_args.modality

def main(analysis_args:Analysis,
         abcd_sync:Path,
         code_dir:Path,
         nperms:int,
         max_perms_per_cpu:int,
         outdir:Path,
         cluster:str,
         mem_free:float,
         run_time:str):
    """Main function"""
    job_dir = outdir /  f'nperms-{nperms}_{analysis_args.recap()}'
    tmp_dir = job_dir / 'tmp'
    tmp_dir.mkdir(parents=True, exist_ok=True)
    log_dir = job_dir / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir/'cluster_submission.log'
    logging.basicConfig(level=logging.DEBUG, filename=log_file, format='%(asctime)s - %(levelname)s:%(message)s')
    print(f"Created log file: {log_file}")

    tabulated_data_dir = abcd_sync / analysis_args.data_release / 'tabulated/released'
    pihat_path = analysis_args.get_genetics(abcd_sync)

    imaging_dir = set_imaging_dir(abcd_sync, analysis_args)

    my_queue = Queue(cluster=cluster, mem_free=mem_free, run_time=run_time)

    n_jobs = nperms//max_perms_per_cpu
    job_ids = []
    for k, perms_list in enumerate(chunks(range(nperms), max_perms_per_cpu), start=1):
        job_nperms = len(perms_list)
        logging.info(f"Starting {job_nperms} permutations ({k}/{n_jobs}).")
        FEMA_wrapper = Matlab(
            'FEMA_wrapper',
            analysis_args.fstem_imaging, analysis_args.design_matrix, tmp_dir/str(k),
            tabulated_data_dir, imaging_dir, analysis_args.datatype,
            **{
                'ico': analysis_args.ico,
                'ranknorm': analysis_args.ranknorm,
                'contrasts': analysis_args.contrasts,
                'RandomEffects': analysis_args.random_effects,
                'mediation': analysis_args.mediation,
                'PermType': analysis_args.permutation_type,
                'colsinterest': analysis_args.colsinterest,
                'tfce': analysis_args.tfce,
                'pihat_file': pihat_path,
                'nperms': job_nperms,
            }
        ).create_command()
        logging.debug(f"matlab_command={FEMA_wrapper}")
        job_wrapper = Job("FEMA_fit", code_dir, FEMA_wrapper)
        job_ids.append(my_queue.submit_job(job_wrapper, dependencies=[], log_dir=log_dir))

    # Gather jobs
    #TODO: add condition if doing MOSTest, TFCE and add arguments
    FEMA_gather = Matlab(
        'FEMA_cluster_gather',
        job_dir,tmp_dir,nperms,**{
            'mostest_reg':'auto'
        }
    ).create_command()
    logging.debug(f"FEMA_gather={FEMA_gather}")
    job_gather = Job("cluster_gather", code_dir, FEMA_gather)
    gather_job_id = my_queue.submit_job(job_gather, dependencies=job_ids, log_dir=log_dir)

    #TODO: if cleanup & concat file exists, remove tmp dir
    logging.info(f"To remove all logs, run:\nrm {log_dir}/*")
    return job_ids, gather_job_id

def parse_args():
    """Argument parser and checker"""
    parser = argparse.ArgumentParser(
        description=(
            'ABCD FEMA Analyses - Cluster submission script.\n'
            'You need to have relevant data stored in a `abcd-sync` directory, access to a job '
            'submission system like SGE or SLURM, a working MATLAB license.'
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    abcd_config_help = "This can be skipped if ~/abcdConfig.json exists with all required paths."

    parse_analysis = parser.add_argument_group('Analysis parameters', 'arguments relevant to the analysis parameters')
    parse_analysis.add_argument('--design_matrix',
        '-dm',
        type=Path,
        help='Path to the design matrix',
        required=True)
    parse_analysis.add_argument('--datatype','-dt',
        type=str,
        help='Choose data type',
        choices=['voxel', 'vertex', 'external', 'corrmat'],
        required=True)
    parse_analysis.add_argument('--modality','-mod',
        help='Choose modality', #TODO: add choices based on abcd-sync options, e.g. for vertex: smri, dmri, tfmri ...
        required=True)
    parse_analysis.add_argument('--fstem_imaging','-img',
        help='Choose imaging type', #TODO: add choices based on abcd-sync options
        required=True)
    parse_analysis.add_argument('--colsinterest','-cols',
        type=str,
        default=None,
        help="Choose which column numbers (in design matrix) to keep for statistical analyses. NB: skip 1-4; so if interested at col 5, input 5-4=[1]") #TODO check if passed only one number vs. several in brackets
    parse_analysis.add_argument('--contrasts',
        default=[])
    parse_analysis.add_argument('--ranknorm',
        default=0,
        help='Whether to rank normalize the imaging data')
    parse_analysis.add_argument('--random_effects','-re',
        default="{'F','S','E'}",
        help='What random effects to include')
    parse_analysis.add_argument('--tfce',
        default=1,
        help='Whether to run TFCE analysis')
    parse_analysis.add_argument('--mediation',
        default=0,
        help='Whether to run mediation analysis')
    parse_analysis.add_argument('--permutation_type','-pt',
        choices=['wildbootstrap', 'wildbootstrap-nn'],
        default='wildbootstrap',
        help=(
            "Which resampling method to use:  "
            "'wildbootstrap' - residual boostrap --> creates null distribution by randomly flipping the sign of each observation; "
            "'wildbootstrap-nn' - non-null boostrap --> estimates distribution around effect of interest using sign flipping (used for sobel test)."
        ))
    parse_analysis.add_argument('--ico',
        choices=[0,1,2,3,4,5,6,7],
        default=5,
        help='ico-number for vertexwise analyses (0-based)')
    parse_analysis.add_argument('--smoothing_factor','-sf',
        default=256)
    parse_analysis.add_argument('--data_release','-v',
        type=str,
        choices=['3.0', '4.0'],
        help='ABCD Data release to use',
        default='4.0')

    parse_dirs = parser.add_argument_group('Directories', 'paths to required for inputs/outputs')
    parse_dirs.add_argument('--outdir','-o',
        type=Path,
        help='Path to the analysis output directory',
        required=True)
    parse_dirs.add_argument('--abcd_sync',
        type=Path,
        help=f'Path to local abcd-sync. {abcd_config_help}',
        required=False)
    parse_dirs.add_argument('--code_dir',
        type=Path,
        help=f'Path to directory containing the relevant code repositories. {abcd_config_help}',
        required=False)
    
    parse_submission = parser.add_argument_group('Cluster info', 'arguments relevant to the cluster submission')
    parse_submission.add_argument('--nperms','-n',
        type=int,
        help='Number of permutations to perform',
        default=10000)
    parse_submission.add_argument('--max_perms_per_cpu','-mp',
        type=int,
        help='Number of permutations to perform per job request',
        default=50) #TODO add max_cpu as an alternative
    parse_submission.add_argument('--mem_free',
        type=float,
        help='Memory (in GB) to request for each job - should depend on type of data analysed.',
        default=20) #TODO precompute a few examples
    parse_submission.add_argument('--run_time',
        type=str,
        help='Run time to request for each job - should depend on type of data analysed (HH:MM:SS).',
        default='02:00:00') #TODO precompute a few examples
    parse_submission.add_argument('--cluster','-c',
        type=str,
        choices=['SLURM', 'SGE', 'LOCAL'],
        help='Which scheduler is installed on your cluster',
        default='SGE')

    args = parser.parse_args()

    if args.abcd_sync is None or args.code_dir is None:
        import json
        abcd_config_path = Path().home() / 'abcdConfig.json'
        try:
            with open(abcd_config_path, 'r') as abcd_config_file:
                abcd_config = json.load(abcd_config_file)
        except FileNotFoundError as err:
            print('Error: ~/abcdConfig.json must exist if not providing `abcd_sync` or `code_dir` arguments')
            raise err
        if args.abcd_sync is None:
            args.abcd_sync = Path(abcd_config['data']['abcd_sync'])
        if args.code_dir is None:
            args.code_dir  = Path(abcd_config['code']['cmig_tools'])
    check_repos(args.code_dir)
    analysis_args = Analysis(
        design_matrix=args.design_matrix,
        datatype=args.datatype,
        modality=args.modality,
        fstem_imaging=args.fstem_imaging,
        colsinterest=args.colsinterest,
        contrasts=args.contrasts,
        ranknorm=args.ranknorm,
        random_effects=args.random_effects,
        tfce=args.tfce,
        mediation=args.mediation,
        permutation_type=args.permutation_type,
        ico=args.ico,
        smoothing_factor=args.smoothing_factor,
        data_release=args.data_release
    )
    return (
        analysis_args,
        args.abcd_sync,
        args.code_dir,
        args.nperms,
        args.max_perms_per_cpu,
        args.outdir.absolute(),
        args.cluster,
        args.mem_free,
        args.run_time
    )


if __name__ == '__main__':
    main(*parse_args())
