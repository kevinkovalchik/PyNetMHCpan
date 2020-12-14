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
from configparser import ConfigParser
from uuid import uuid4


ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
config_file = str(Path(ROOT_DIR)/'pynetmhcpan.config')
common_aa = "ARNDCQEGHILKMFPSTWYV"

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


def remove_modifications(peptides: Union[List[str], str]):
    if isinstance(peptides, str):
        return ''.join(re.findall('[a-zA-Z]+', peptides))
    unmodified_peps = []
    for pep in peptides:
        pep = ''.join(re.findall('[a-zA-Z]+', pep))
        unmodified_peps.append(pep)
    return unmodified_peps


def remove_previous_and_next_aa(peptides: Union[List[str], str]):
    return_one = False
    if isinstance(peptides, str):
        peptides = [peptides]
        return_one = True
    for i in range(len(peptides)):
        if peptides[i][1] == '.':
            peptides[i] = peptides[i][2:]
        if peptides[i][-2] == '.':
            peptides[i] = peptides[i][:-2]
    if return_one:
        return peptides[0]
    return peptides


def replace_uncommon_aas(peptide):
    pep = peptide
    for aa in peptide:
        if aa not in common_aa:
            pep = pep.replace(aa, 'X')
    return pep


def create_netmhcpan_peptide_index(peptide_list):
    netmhcpan_peps = {}
    for i in range(len(peptide_list)):
        if len(peptide_list[i]) < 1:
            continue
        netmhc_pep = replace_uncommon_aas(peptide_list[i])
        netmhcpan_peps[peptide_list[i]] = netmhc_pep
    return netmhcpan_peps


class Helper:
    """
    example usage:
    cl_tools.make_binding_prediction_jobs()
    cl_tools.run_jubs()
    cl_tools.aggregate_netmhcpan_results()
    cl_tools.clear_jobs()
    """
    def __init__(self,
                 peptides: List[str] = None,
                 mhc_class: str = 'I',
                 alleles: List[str] = ('HLA-A03:02', 'HLA-A02:02'),
                 n_threads: int = 0,
                 tmp_dir: str = '/tmp',
                 output_dir: str = None):
        """
        Helper class to run NetMHCpan or NetMHCIIpan on multiple CPUs from Python.

        :param peptides: A list of peptides.
        :param mhc_class: The MHC class of the peptides (i.e. I or II)
        :param alleles: A list of NetMHCpan- or NetMHCIIpan-recognized alleles.
        :param min_length: The minimum peptide length which will be analyzed. Minimum possible values are
        8 for class I and 9 for class II.
        :param max_length: The maximum peptide length to analyze.
        :param netmhcpan_version: What version of NetMHCpan is being used. Option are 4.0 and 4.1.
        :param netmhcpan_executable: Location of the NetMHCpan execution script. The default value assumes it is
        in your PATH.
        :param netmhc2pan_executable: Location of the NetMHCIIpan execution script. The default value assumes it is
        in your PATH.
        :param n_threads: The number of threads to use
        :param output_dir: The output or temp directory. If you are going to use the analysis results programmatically
        and don't need the output files you can set this to something like /temp/netmhcpan
        """

        config = ConfigParser()
        config.read(config_file)

        self.NETMHCPAN = config['PATHS']['NetMHCpan']
        self.NETMHCIIPAN = config['PATHS']['NetMHCIIpan']

        if isinstance(alleles, str):
            if ',' in alleles:
                alleles = alleles.split(',')
            elif ' ' in alleles:
                alleles = alleles.split(' ')
            else:
                alleles = [alleles]
        self.alleles = alleles
        self.mhc_class = mhc_class
        self.min_length = 8 if self.mhc_class == 'I' else 9
        self.peptides = peptides
        self.netmhcpan_peptides = dict()
        self.predictions = pd.DataFrame(
            columns=['Peptide', 'Allele', 'Rank', 'Binder']
        )
        self.wd = Path(output_dir) if output_dir else Path(os.getcwd())
        self.temp_dir = Path(tmp_dir) / 'PyNetMHCpan'
        if not self.wd.exists():
            self.wd.mkdir(parents=True)
        if not self.temp_dir.exists():
            self.temp_dir.mkdir(parents=True)
        self.predictions_made = False
        self.not_enough_peptides = []
        if n_threads < 1 or n_threads > os.cpu_count():
            self.n_threads = os.cpu_count()
        else:
            self.n_threads = n_threads
        self.jobs = []

        if (self.mhc_class == 'I' and not Path(self.NETMHCPAN).is_file()) or \
                (self.mhc_class == 'II' and not Path(self.NETMHCIIPAN).is_file()):
            print("It looks like there is a problem with the configuration. Please review the "
                  "following and update as necessary:")
            self.update_config()

    def update_config(self):
        with open(config_file, 'r') as f:
            config = f.readlines()

        print('The current config settings are as follows:')
        print()
        for line in config:
            print(f' {line}', end='')
        print()
        netmhcpan = input('To update the path to the NetMHCpan script type 1 and press enter '
                          '(any other key to skip): ')

        if netmhcpan == '1':
            new_netmhcpan = input('Type the new path to the NetMHCpan execution script and press enter: ')
        else:
            netmhcpan = None
            new_netmhcpan = None
        netmhc2pan = input('To update the path to the NetMHCIIpan script type 2 and press enter '
                           '(any other key to skip): ')
        if netmhc2pan == '2':
            new_netmhc2pan = input('Type the new path to the NetMHCIIpan execution script and press enter: ')
        else:
            netmhc2pan = None
            new_netmhc2pan = None

        if netmhcpan or netmhc2pan:
            config = ConfigParser()
            config.read(config_file)

            with open(config_file, 'w') as f:
                f.write('[PATHS]\n')
                if netmhcpan:
                    f.write(f'NetMHCpan = {new_netmhcpan}\n')
                else:
                    f.write(f'NetMHCpan = {config["PATHS"]["NetMHCpan"]}\n')
                if netmhc2pan:
                    f.write(f'NetMHCIIpan = {new_netmhc2pan}\n')
                else:
                    f.write(f'NetMHCIIpan = {config["PATHS"]["NetMHCIIpan"]}\n')

            config = ConfigParser()
            config.read(config_file)
            # update the paths
            self.NETMHCPAN = config['PATHS']['NetMHCpan']
            self.NETMHCIIPAN = config['PATHS']['NetMHCIIpan']

    def add_peptides(self, peptides: List[str]):
        if not self.peptides:
            self.peptides = []
        peptides = remove_previous_and_next_aa(peptides)
        peptides = remove_modifications(peptides)
        use_peptides = [x for x in peptides if len(x) >= self.min_length]
        self.peptides += use_peptides

        if len(use_peptides) < len(peptides):
            print(f'{len(peptides)-len(use_peptides)} peptides were outside the length restrictions and '
                  f'were removed from the peptide list.')

        self.netmhcpan_peptides = create_netmhcpan_peptide_index(self.peptides)

    def _make_binding_prediction_jobs(self):
        if not self.peptides:
            print("ERROR: You need to add some peptides first!")
            return
        self.jobs = []

        # split peptide list into chunks
        peptides = np.array(list(self.netmhcpan_peptides.values()))
        np.random.shuffle(peptides)  # we need to shuffle them so we don't end up with files filled with peptide lengths that take a LONG time to compute (this actually is a very significant speed up)

        if len(peptides) > 100:
            chunks = array_split(peptides, self.n_threads)
        else:
            chunks = [peptides]
        job_number = 1

        for chunk in chunks:
            if len(chunk) < 1:
                continue
            fname = Path(self.temp_dir, f'peplist_{job_number}.csv')
            # save the new peptide list, this will be given to netMHCpan
            chunk.tofile(str(fname), '\n', '%s')
            # run netMHCpan
            if self.mhc_class == 'I':
                command = f'{self.NETMHCPAN} -p -f {fname} -a {",".join(self.alleles)}'.split(' ')
            else:
                command = f'{self.NETMHCIIPAN} -inptype 1 -f {fname} -a {",".join(self.alleles)}'.split(' ')
            job = Job(command=command,
                      working_directory=self.temp_dir)
            self.jobs.append(job)
            job_number += 1

    def _run_jobs(self):
        self.jobs = _run_multiple_processes(self.jobs, n_processes=self.n_threads)

    def _clear_jobs(self):
        self.jobs = []

    def _aggregate_netmhcpan_results(self):
        for job in self.jobs:
            self._parse_netmhc_output(job.stdout.decode())

        self.predictions.to_csv(str(Path(self.temp_dir) / f'netMHC'
                                                          f'{"II" if self.mhc_class == "II" else ""}'
                                                          f'pan_predictions.csv'))

    def _parse_netmhc_output(self, stdout: str):
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

    def make_predictions(self):
        self.temp_dir = self.temp_dir / str(uuid4())
        self.temp_dir.mkdir(parents=True)
        self._make_binding_prediction_jobs()
        self._run_jobs()
        self._aggregate_netmhcpan_results()
        self._clear_jobs()

    def annotate_file(self, filename: str,
                      peptide_column: str = 'Peptide',
                      delimiter: str = '\t'):

        # clear the peptide list
        self.peptides = []

        with open(filename, 'r') as f:
            header = f.readline().strip().split(delimiter)
            content = [x.strip().split(delimiter) for x in f.readlines()]
        pep_index = header.index(peptide_column)
        peptides = [x[pep_index] for x in content]
        self.add_peptides(peptides)
        self.make_predictions()
        binding = {p: {} for p in self.predictions['Peptide'].unique()}
        for i in self.predictions.index:
            pep, allele, rank, binder = self.predictions.loc[i, :]
            binding[pep][allele] = rank

        alleles = self.predictions['Allele'].unique()
        new_content = []

        for allele in alleles:
            header.insert(pep_index-1, f'{allele}_rank')
        for line in content:
            pep = replace_uncommon_aas(remove_modifications(remove_previous_and_next_aa(line[pep_index])))
            # pep = self.netmhcpan_peptides[pep]
            if len(pep) < self.min_length:
                continue
            ranks = []
            for allele in alleles:
                rank = float(binding[pep][allele])
                #line.insert(pep_index-1, str(np.log(rank)))
                line.insert(pep_index - 1, str(rank))
                ranks.append(rank)
            new_content.append(line)

        f_out = self.wd / (Path(filename).stem + '_annotated.tsv')
        with open(f_out, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for line in new_content:
                f.write('\t'.join(line) + '\n')



