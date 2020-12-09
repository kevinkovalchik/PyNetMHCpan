import subprocess
import pandas as pd
import os
import re
from numpy import array_split
import numpy as np
from pathlib import Path
from typing import List, Union
from multiprocessing import Pool
from datetime import datetime


class Job:
    def __init__(self,
                 command: Union[str, List[str]],
                 working_directory: Union[str, Path, None],
                 sample=None):

        self.command = command
        self.working_directory = working_directory
        self.returncode = None
        self.time_start = str(datetime.now()).replace(' ', '')
        self.time_end = ''
        self.stdout = ''
        self.stderr = ''
        self.sample = sample

    def run(self):
        if self.working_directory is not None:
            os.chdir(self.working_directory)

        command = self.command.split(' ') if isinstance(self.command, str) else self.command
        p = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        self.stdout, self.stderr = p.communicate()
        self.time_end = str(datetime.now()).replace(' ', '')
        self.returncode = p.returncode


def run(job: Job):
    job.run()
    return job


def _run_multiple_processes(jobs: List[Job], n_processes: int):
    pool = Pool(n_processes)
    returns = pool.map(run, jobs)
    pool.close()
    return returns


def remove_modifications(peptide_list, verbose=False):
    unmodified_peps = []
    if verbose:
        print('Removing peptide modifications')
    for pep in peptide_list:
        pep = ''.join(re.findall('[a-zA-z]+', pep))
        unmodified_peps.append(pep)
    return unmodified_peps


class MhcToolHelper:
    """
    example usage:
    cl_tools.make_binding_prediction_jobs()
    cl_tools.run_jubs()
    cl_tools.aggregate_netmhcpan_results()
    cl_tools.clear_jobs()
    """
    def __init__(self,
                 peptides: List[str],
                 mhc_class: str = 'I',
                 alleles: List[str] = ('HLA-A03:02', 'HLA-A02:02'),
                 min_length: int = 8,
                 max_length: int = 12,
                 netmhcpan_version: float = 4.1,
                 netmhcpan_executable: str = "netMHCpan",
                 netmhc2pan_executable: str = "netMHCIIpan",
                 n_threads: int = 0,
                 output_dir: str = None):

        if mhc_class == 'I' and min_length < 8:
            raise ValueError('Class I peptides must be 8 mers and longer for NetMHCpan')
        if mhc_class == 'II' and min_length < 9:
            raise ValueError('Class II peptides must be 9 mers and longer for NetMHCIIpan')

        self.NETMHCPAN = netmhcpan_executable
        self.NETMHCIIPAN = netmhc2pan_executable
        self.NETMHCPAN_VERSION = netmhcpan_version

        if isinstance(alleles, str):
            if ',' in alleles:
                alleles = alleles.split(',')
            elif ' ' in alleles:
                alleles = alleles.split(' ')
            else:
                alleles = [alleles]
        self.alleles = alleles
        self.mhc_class = mhc_class
        self.min_length = min_length
        self.max_length = max_length
        self.peptides = peptides
        self.predictions = pd.DataFrame(
            columns=['Peptide', 'Allele', 'Rank', 'Binder']
        )
        self.wd = Path(output_dir) if output_dir else os.getcwd()
        if not self.wd.exists():
            self.wd.mkdir(parents=True)
        self.predictions_made = False
        self.not_enough_peptides = []
        if n_threads < 1 or n_threads > os.cpu_count():
            self.n_threads = os.cpu_count()
        else:
            self.n_threads = n_threads
        self.jobs = []

        with open(str(self.wd / 'all_peptides.tsv'), 'w') as f:
            for pep in peptides:
                f.write(pep + '\n')

    def make_binding_prediction_jobs(self):
        # split peptide list into chunks
        peptides = np.array(remove_modifications(self.peptides))
        lengths = np.vectorize(len)(peptides)
        peptides = peptides[(lengths >= self.min_length) & (lengths <= self.max_length)]
        np.random.shuffle(peptides)  # we need to shuffle them so we don't end up with files filled with peptide lengths that take a LONG time to compute (this actually is a very significant speed up)

        if len(peptides) > 100:
            chunks = array_split(peptides, self.n_threads)
        else:
            chunks = [peptides]
        job_number = 1
        results = []

        for chunk in chunks:
            if len(chunk) < 1:
                continue
            fname = Path(self.wd, f'peplist_{job_number}.csv')
            # save the new peptide list, this will be given to netMHCpan
            chunk.tofile(str(fname), '\n', '%s')
            # run netMHCpan
            if self.mhc_class == 'I':
                command = f'{self.NETMHCPAN} -p -f {fname} -a {",".join(self.alleles)}'.split(' ')
            else:
                command = f'{self.NETMHCIIPAN} -inptype 1 -f {fname} -a {",".join(self.alleles)}'.split(' ')
            job = Job(command=command,
                      working_directory=self.wd)
            self.jobs.append(job)
            job_number += 1

    def run_jobs(self):
        self.jobs = _run_multiple_processes(self.jobs, n_processes=self.n_threads)

    def clear_jobs(self):
        self.jobs = []

    def aggregate_netmhcpan_results(self):
        for job in self.jobs:
            self.parse_netmhc_output(job.stdout.decode())

        self.predictions.to_csv(str(Path(self.wd) / f'netMHC'
                                                          f'{"II" if self.mhc_class == "II" else ""}'
                                                          f'pan_predictions.csv'))

    def parse_netmhc_output(self, stdout: str):
        rows = []
        lines = stdout.split('\n')
        if self.mhc_class == 'I':  # works for 4.0 and 4.1, will need to keep an eye on future releases
            allele_idx = 1
            peptide_idx = 2
            rank_idx = 12
        else:  # works for NetMHCIIpan4.0
            allele_idx = 1
            peptide_idx = 2
            rank_idx = 8
        for line in lines:
            line = line.strip()
            line = line.split()
            if not line or line[0] == '#' or not line[0].isnumeric():
                continue
            allele = line[allele_idx].replace('*', '')
            peptide = line[peptide_idx]
            rank = line[rank_idx]
            if self.mhc_class == 'I':
                if float(rank) <= 0.5:
                    binder = 'Strong'
                elif float(rank) <= 2.0:
                    binder = 'Weak'
                else:
                    binder = 'Non-binder'
            else:
                if float(rank) <= 2:
                    binder = 'Strong'
                elif float(rank) <= 10:
                    binder = 'Weak'
                else:
                    binder = 'Non-binder'
            rows.append((peptide, allele, rank, binder))
        self.predictions = self.predictions.append(
            pd.DataFrame(columns=['Peptide', 'Allele', 'Rank', 'Binder'], data=rows),
            ignore_index=True
        )
        if len(rows) == 0:
            print(stdout)

    def get_predictions(self):
        self.make_binding_prediction_jobs()
        self.run_jobs()
        self.aggregate_netmhcpan_results()
        self.clear_jobs()
