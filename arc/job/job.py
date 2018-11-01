#!/usr/bin/env python
# encoding: utf-8

import os
import time
import csv
import logging

from arc.settings import arc_path, software_server, servers, submit_filename, output_filename
from arc.job.submit import submit_sctipts
from arc.job.input import input_files
from arc.job.ssh import SSH_Client

##################################################################


class Job(object):
    """
    ARC Job class. The attributes are:

    ================ =================== ===============================================================================
    Attribute        Type                Description
    ================ =================== ===============================================================================
    `project`         ``str``            The project's name. Used for naming the directory.
    `species_name`    ``str``            The species/TS name. Used for naming the directory.
    `charge`          ``int``            The species net charge. Default is 0
    `multiplicity`    ``int``            The species multiplicity.
    `spin`            ``int``            The spin. automatically derived from the multiplicity
    `xyz`             ``str``            The xyz geometry. Used for the calculation
    `conformer`       ``int``            Conformer number if optimizing conformers
    `is_ts`           ``bool``           Whether this species represents a transition structure
    `level_of_theory` ``str``            Level of theory, e.g. 'CBS-QB3', 'CCSD(T)-F12a/aug-cc-pVTZ',
                                           'B3LYP/6-311++G(3df,3pd)'...
    `job_type`         ``str``           The job's type
    `scan`             ``str``           A string representing atom labels for the dihedral scan (e.g., '2, 1, 3, 5')
    `software`         ``str``           The electronic structure software to be used
    `memory`           ``int``           The allocated memory (1000 mb by default)
    `method`           ``str``           The calculation method (e.g., 'B3LYP', 'CCSD(T)', 'CBS-QB3'...)
    `basis_set`        ``str``           The basis set (e.g., '6-311++G(d,p)', 'aug-cc-pVTZ'...)
    `fine`             ``bool``          Whether to use fine geometry optimization parameters
    `shift`            ``str``           A string representation alpha- and beta-spin orbitals shifts (molpro only)
    `comments`         ``str``           Job comments (archived, not used)
    `date`             ``str``           The date this job was initiated. Determined automatically
    `run_time`         ``float``         Runtime is seconds. Determined automatically
    `job_status`       ``list``          The job's server and ESS statuses. Determined automatically
    `job_name`         ``str``           Job's name (e.g., 'a103'). Determined automatically
    `job_id`           ``int``           The job's ID determined by the server.
    `local_path`       ``str``           Local path to job's folder. Determined automatically
    `remote_path`      ``str``           Remote path to job's folder. Determined automatically
    `submit`           ``str``           The submit script. Created automatically
    `input`            ``str``           The input file. Created automatically
    `server`           ``str``           Server's name. Determined automatically
    ================ =================== ===============================================================================
    """
    def __init__(self, project, species_name, xyz, job_type, level_of_theory, multiplicity, charge=0, conformer=-1,
                 fine=False, shift='', software=None, is_ts=False, scan='', memory=1000, comments=''):
        self.project = project
        self.species_name = species_name
        self.job_num = -1
        self.charge = charge
        self.multiplicity = multiplicity
        self.spin = self.multiplicity - 1
        self.xyz = xyz
        self.conformer = conformer
        self.is_ts = is_ts
        job_types = ['opt', 'ts', 'freq', 'optfreq', 'tsfreq', 'sp', 'composite', 'scan', 'irc', 'gsm']
        if job_type not in job_types:
            raise ValueError("Job type {0} not understood. Must be on of the following: {1}".format(
                job_type, job_types))
        self.job_type = job_type

        # determine level of theory and software to use:
        self.level_of_theory = level_of_theory.lower()
        self.method, self.basis_set = '', ''
        if job_type == 'composite':
            self.method, self.basis_set = self.level_of_theory, ''
        else:
            self.method, self.basis_set = self.level_of_theory.split('/')
        if software is not None:
            self.software = software.lower()
        else:
            if job_type == 'composite':
                self.software = 'gaussian03'
            elif job_type in ['opt', 'ts', 'optfreq', 'tsfreq', 'sp']:
                if 'ccs' in self.method or 'cis' in self.method or 'pv' in self.basis_set:
                    self.software = 'molpro2012'
                elif 'b3lyp' in self.method:
                    self.software = 'gaussian03'
                elif 'wB97' in self.method or 'm06-2x' in self.method:
                    self.software = 'qchem'
            elif job_type == 'scan':
                self.software = 'gaussian03'
            elif job_type == 'gsm':
                self.software = 'gaussian03'
            elif job_type == 'irc':
                self.software = 'gaussian03'

        self.scan = scan
        self.memory = memory
        self.fine = fine
        self.shift = shift
        self.date = time.asctime()
        self.run_time = 0
        self.job_status = ['initializing', 'initializing']
        self.job_id = 0
        self.comments = comments

        self._set_job_number()
        self.job_name = 'a' + str(self.job_num)
        conformer_folder = '' if self.conformer < 0 else os.path.join('conformers',str(self.conformer))
        self.local_path = os.path.join(arc_path, 'Projects', self.project,
                                       self.species_name, conformer_folder, self.job_name)
        self.remote_path = os.path.join('runs', 'ARC_Projects', self.project,
                                       self.species_name, conformer_folder, self.job_name)
        self.submit = ''
        self.input = ''

        self.server = software_server[self.software]
        self._write_initiated_job_to_csv_file()

    def _set_job_number(self):
        """
        Used only as the entry number in csv archiving
        """
        csv_path = os.path.join(arc_path, 'initiated_jobs.csv')
        if not os.path.isfile(csv_path):
            # check file, make index file and write headers if file doesn't exists
            with open(csv_path, 'wb') as f:
                writer = csv.writer(f, dialect='excel')
                row = ['job_num', 'project', 'species_name', 'conformer', 'is_ts', 'charge', 'multiplicity', 'job_type',
                       'job_name', 'job_id', 'server', 'software', 'memory', 'method', 'basis_set', 'date', 'comments']
                writer.writerow(row)
        with open(csv_path, 'rb') as f:
            reader = csv.reader(f, dialect='excel')
            job_num = 0
            for row in reader:
                job_num += 1
            self.job_num = job_num

    def _write_initiated_job_to_csv_file(self):
        """
        Write an initiated ARCJob into the initiated_jobs.csv file.
        """
        csv_path = os.path.join(arc_path, 'initiated_jobs.csv')
        if self.conformer < 0:  # this is not a conformer search job
            conformer = '-'
        else: conformer = str(self.conformer)
        with open(csv_path, 'ab') as f:
            writer = csv.writer(f, dialect='excel')
            row = [self.job_num, self.project, self.species_name, conformer, self.is_ts, self.charge,
                   self.multiplicity, self.job_type, self.job_name, self.job_id, self.server, self.software,
                   self.memory, self.method, self.basis_set, self.date, self.comments]
            writer.writerow(row)

    def write_completed_job_to_csv_file(self):
        """
        Write a completed ARCJob into the completed_jobs.csv file.
        """
        csv_path = os.path.join(arc_path, 'completed_jobs.csv')
        if not os.path.isfile(csv_path):
            # check file, make index file and write headers if file doesn't exists
            with open(csv_path, 'wb') as f:
                writer = csv.writer(f, dialect='excel')
                row = ['job_num', 'project', 'species_name', 'conformer', 'is_ts', 'charge', 'multiplicity', 'job_type',
                       'job_name', 'job_id', 'server', 'software', 'memory', 'method', 'basis_set', 'date', 'run_time',
                       'job_status_(server)', 'job_status_(ESS)', 'comments']
                writer.writerow(row)
        csv_path = os.path.join(arc_path, 'completed_jobs.csv')
        if self.conformer < 0:  # this is not a conformer search job
            conformer = '-'
        else: conformer = str(self.conformer)
        with open(csv_path, 'ab') as f:
            writer = csv.writer(f, dialect='excel')
            row = [self.job_num, self.project, self.species_name, conformer, self.is_ts, self.charge,
                   self.multiplicity, self.job_type, self.job_name, self.job_id, self.server, self.software,
                   self.memory, self.method, self.basis_set, self.date, self.run_time, self.job_status[0],
                   self.job_status[1], self.comments]
            writer.writerow(row)

    def write_submit_script(self):
        un = servers[self.server]['un']  # user name
        self.submit = submit_sctipts[self.software].format(self.job_name, un)
        if not os.path.exists(self.local_path):
            os.makedirs(self.local_path)
        with open(os.path.join(self.local_path, submit_filename[self.server]), 'wb') as f:
            f.write(self.submit)
        self._upload_submit_file()

    def write_input_file(self):
        """
        Write a software-specific job-specific input file.
        Saves the file locally and also uploads it to the server.
        """

        self.input = input_files[self.software]

        slash = ''
        if self.software == 'gaussian03' and not self.job_type == 'composite':
            slash = '/'

        if self.multiplicity > 1:
            restricted = 'u'
        else:
            restricted = ''

        job_type_1, job_type_2, fine = '', '', ''
        if self.job_type == 'opt':
            if self.software == 'gaussian03':
                job_type_1 = 'opt=calcfc'
            elif self.software == 'qchem':
                job_type_1 = 'opt'
                if self.fine:
                    fine = '\n   GEOM_OPT_TOL_GRADIENT 15\n   GEOM_OPT_TOL_DISPLACEMENT 60\n   GEOM_OPT_TOL_ENERGY 5'
            elif 'molpro' in self.software:
                job_type_1 = "\noptg,savexyz='geometry.xyz'"

        if self.job_type == 'ts':
            if self.software == 'gaussian03':
                job_type_1 = 'opt=(calcfc,ts,noeigen)'
            elif self.software == 'qchem':
                job_type_1 = 'ts'
                if self.fine:
                    fine = '\n   GEOM_OPT_TOL_GRADIENT 15\n   GEOM_OPT_TOL_DISPLACEMENT 60\n   GEOM_OPT_TOL_ENERGY 5'
            elif 'molpro' in self.software:
                job_type_1 = "\noptg,root=2,method=qsd,readhess,savexyz='geometry.xyz'"

        elif self.job_type == 'freq':
            if self.software == 'gaussian03':
                job_type_2 = 'freq iop(7/33=1)'
            elif self.software == 'qchem':
                job_type_1 = 'freq'
            elif 'molpro' in self.software:
                job_type_1 = '\n{frequencies;\nthermo;\nprint,HESSIAN,thermo;}'

        elif self.job_type == 'optfreq':
            if self.software == 'gaussian03':
                job_type_1 = 'opt=calcfc'
                job_type_2 = 'freq iop(7/33=1)'
            elif self.software == 'qchem':
                self.input += """@@@

$molecule
   read
$end

$rem
   JOBTYPE       {job_type_2}
   METHOD        {method}
   UNRESTRICTED  {restricted}
   BASIS         {basis}
   SCF_GUESS     read
$end

"""
                job_type_1 = 'opt'
                job_type_2 = 'freq'
                if self.fine:
                    fine = '\n   GEOM_OPT_TOL_GRADIENT 15\n   GEOM_OPT_TOL_DISPLACEMENT 60\n   GEOM_OPT_TOL_ENERGY 5'
            elif 'molpro' in self.software:
                job_type_1 = "\noptg,savexyz='geometry.xyz"
                job_type_2 = '\n{frequencies;\nthermo;\nprint,HESSIAN,thermo;}'

        elif self.job_type == 'tsfreq':
            if self.software == 'gaussian03':
                job_type_1 = 'opt=(calcfc,ts,noeigen)'
                job_type_2 = 'freq iop(7/33=1)'
            elif self.software == 'qchem':
                self.input += """@@@

$molecule
   read
$end

$rem
   JOBTYPE       {job_type_2}
   METHOD        {method}
   UNRESTRICTED  {restricted}
   BASIS         {basis}
   SCF_GUESS     read
$end

"""
                job_type_1 = 'ts'
                job_type_2 = 'freq'
                if self.fine:
                    fine = '\n   GEOM_OPT_TOL_GRADIENT 15\n   GEOM_OPT_TOL_DISPLACEMENT 60\n   GEOM_OPT_TOL_ENERGY 5'
            elif 'molpro' in self.software:
                job_type_1 = "\noptg,root=2,method=qsd,readhess,savexyz='geometry.xyz'"
                job_type_2 = '\n{frequencies;\nthermo;\nprint,HESSIAN,thermo;}'

        if self.job_type == 'sp':
            if self.software == 'gaussian03':
                pass
            elif self.software == 'qchem':
                job_type_1 = 'sp'
            elif 'molpro' in self.software:
                pass

        if self.job_type == 'composite':
            if self.software == 'gaussian03':
                if self.is_ts:
                    job_type_1 = 'opt=(ts,noeigentest,calcfc)'
            else:
                raise ValueError("Currently composite methods are only available in gaussian03")

        if self.job_type == 'scan':
            if self.software == 'gaussian03':
                job_type_1 = 'opt=modredundant'
                self.scan = '\nD ' + self.scan + ' S 36 10.000000'
            else:
                raise ValueError("Currently rotor scan is only available in gaussian03")

        if self.job_type == 'irc':  # TODO
            pass

        if self.job_type == 'gsm':  # TODO
            pass

        self.input = self.input.format(memory=self.memory, method=self.method, slash=slash, basis=self.basis_set,
                                       charge=self.charge, multiplicity=self.multiplicity, spin=self.spin, xyz=self.xyz,
                                       job_type_1=job_type_1, job_type_2=job_type_2, scan=self.scan,
                                       restricted=restricted, fine=fine, shift=self.shift)

        if not os.path.exists(self.local_path):
            os.makedirs(self.local_path)
        with open(os.path.join(self.local_path, 'input.in'), 'wb') as f:
            f.write(self.input)
        self._upload_input_file()

    def _upload_submit_file(self):
        ssh = SSH_Client(self.server)
        ssh.send_command_to_server(command='mkdir -p {0}'.format(self.remote_path))
        remote_file_path = os.path.join(self.remote_path, submit_filename[self.server])
        ssh.upload_file(remote_file_path=remote_file_path, file_string=self.submit)

    def _upload_input_file(self):
        ssh = SSH_Client(self.server)
        ssh.send_command_to_server(command='mkdir -p {0}'.format(self.remote_path))
        remote_file_path = os.path.join(self.remote_path, 'input.in')
        ssh.upload_file(remote_file_path=remote_file_path, file_string=self.input)

    def _download_output_file(self):
        ssh = SSH_Client(self.server)
        remote_file_path = os.path.join(self.remote_path, output_filename[self.software])
        local_file_path = os.path.join(self.local_path, 'output.out')
        ssh.download_file(remote_file_path=remote_file_path, local_file_path=local_file_path)

    def run(self):
        logging.info('\nrunning job {0}'.format(self.job_name))
        logging.info('writing submit script...')
        self.write_submit_script()
        logging.info('writing input file...')
        self.write_input_file()
        ssh = SSH_Client(self.server)
        logging.info('submitting job...')
        self.job_status[0], self.job_id = ssh.submit_job(remote_path=self.remote_path)

    def determine_job_status(self):
        if self.job_status[0] == 'errored':
            return
        server_status = self._check_job_server_status()
        ess_status = ''
        if server_status == 'done':
            ess_status = self._check_job_ess_status()
        if server_status == 'running':
            ess_status = 'running'
        self.job_status = [server_status, ess_status]

    def _check_job_server_status(self):
        """
        Possible statuses: `initializing`, `running`, `errored on node xx`, `done`
        """
        ssh = SSH_Client(self.server)
        return ssh.check_job_status(self.job_id)

    def _check_job_ess_status(self):
        """
        Check the status of the job ran by the electronic structure software (ESS)
        Possible statuses: `initializing`, `running`, `errored: {error type / message}`, `unconverged`, `done`
        """
        output_path = os.path.join(self.local_path, 'output.out')
        if os.path.exists(output_path):
            os.remove(output_path)
        self._download_output_file()
        with open(output_path, 'rb') as f:
            lines = f.readlines()
            if self.software == 'gaussian03':
                if 'Normal termination of Gaussian' in lines[-2]:
                    return 'done'
                else:
                    for line in lines[::-1]:
                        if 'Error' in line or 'NtrErr' in line or 'Erroneous' in line or 'malloc' in line\
                                or 'galloc' in line:
                            reason = ''
                            if 'l9999.exe' in line or 'l103.exe' in line:
                                return 'unconverged'
                            if 'l502.exe' in line:
                                return 'unconverged SCF'
                            if 'Erroneous write' in line or 'Write error in NtrExt1' in line:
                                reason = 'Ran out of disk space.'
                            if 'l716.exe' in line:
                                reason = 'Angle in z-matrix outside the allowed range 0 < x < 180.'
                            if 'l301.exe' in line:
                                reason = 'Input Error. Either charge, multiplicity, or basis set was not specified ' \
                                         'correctly. Or, an atom specified  does not match any standard atomic symbol.'
                            if 'NtrErr Called from FileIO' in line:
                                reason = 'Operation on .chk file was specified, but .chk was not found.'
                            if 'l101.exe' in line:
                                reason = 'Input Error. The blank line after the coordinate section is missing, ' \
                                         'or charge/multiplicity was not specified correctly.'
                            if 'l202.exe' in line:
                                reason = 'During the optimization process, either the standard orientation ' \
                                         'or the point group of the molecule has changed.'
                            if 'l401.exe' in line:
                                reason = 'The projection from the old to the new basis set has failed.'
                            if 'malloc failed' in line or 'galloc' in line:
                                reason = 'Memory allocation failed (did you ask for too much?)'
                            if 'A SYNTAX ERROR WAS DETECTED' in line:
                                reason = 'Check .inp carefully for syntax errors in keywords.'
                            return 'errored: {0}; {1}'.format(line, reason)
                    return 'errored: Unknown reason'
            elif self.software == 'qchem':
                for line in lines[-1:-15]:
                    if 'opt' in self.job_type or 'ts' in self.job_type:
                        for line in lines[::-1]:
                            if 'MAXIMUM OPTIMIZATION CYCLES REACHED' in line:
                                return 'unconverged'
                            if 'OPTIMIZATION CONVERGED' in line:
                                break
                    if 'Thank you very much for using Q-Chem' in line:
                        return 'done'
                for line in lines[::-1]:
                    if 'error' in line:
                        return 'errored: ' + line
                return 'errored: Unknown reason'
            elif 'molpro' in self.software:
                for line in lines[::-1]:
                    if 'Molpro calculation terminated' in line:
                        return 'done'
                    if 'No convergence' in line:
                        return 'unconverged'
                for line in lines[::-1]:
                    if 'the problem occurs' in line:
                        return 'errored: ' + line
                return 'errored: Unknown reason'

    def troubleshoot_server(self):
        # change node? forbid a node on pharos? this method should also delete and resubmit
        pass

# TODO: irc, gsm input files
# TODO: write server troubleshooting method