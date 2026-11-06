# SPDX-License-Identifier: GPL-3.0-or-later
import os
import subprocess
import sys
import tempfile

from .configure import exe


def _amrfinder_env():
    """Point AMRFinderPlus at the database shipped with the conda environment."""
    env = os.environ.copy()
    if 'AMRFINDER_DB' not in env:
        candidate = os.path.join(sys.prefix, 'share', 'amrfinderplus', 'data', 'latest')
        if os.path.exists(candidate):
            env['AMRFINDER_DB'] = candidate
    return env


def run_amrfinder(input, threads):
    with tempfile.TemporaryDirectory(prefix='staphmix_amr_') as tmpdir:
        output_name = os.path.join(tmpdir, 'amrfinder.tsv')
        cmd = [exe['amrfinder'],
               '-n', input,
               '--threads', str(threads),
               '-o', output_name]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True, env=_amrfinder_env())
        if proc.returncode != 0:
            raise RuntimeError('AMRFinderPlus failed:\n{0}'.format(
                (proc.stderr or proc.stdout).strip()))

        data = []
        amr_class = set()
        with open(output_name) as fin:
            fin.readline()
            for line in fin:
                part = line.rstrip('\n').split('\t')
                if len(part) < 19 or part[8] != 'AMR':
                    continue
                data.append(dict(
                    locus_tag='-',
                    contig=part[1],
                    start=int(part[2]),
                    end=int(part[3]),
                    direct=part[4],
                    identity=part[16],
                    coverage=part[15],
                    accession=part[18],
                    gene=part[5],
                    category=part[8] if part[8] == part[9] else '{0}:{1}'.format(part[8], part[9]),
                    subcategory=part[10] if part[10] == part[11] else '{0}:{1}'.format(part[10], part[11]),
                    function=part[6],
                ))
                amr_class.add(part[10])
    return data, len(amr_class), len(data)
