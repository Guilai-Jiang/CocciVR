# SPDX-License-Identifier: GPL-3.0-or-later
"""Unit tests that do not need the reference database or external tools."""

from click.testing import CliRunner

from coccivr import __version__
from coccivr.cli import main
from coccivr.read_seq import fastq_to_fasta
from coccivr.run_score import calculate_amr, calculate_virulence


def test_fastq_to_fasta_keeps_header_marker():
    content = "@read1\nACGTACGT\n+\nIIIIIIII\n"
    assert fastq_to_fasta(content) == ">read1\nACGTACGT\n"


def test_enterotoxin_score():
    hits = [{'gene': 'sea'}, {'gene': 'seb'}, {'gene': 'hlgC'}]
    assert calculate_virulence(hits) == [2, 0]


def test_blood_infection_score():
    hits = [{'gene': 'lukS-PV'}, {'gene': 'tsst-1'}]
    assert calculate_virulence(hits) == [0, 4]


def test_amr_mrsa():
    mec, score, special = calculate_amr([{'gene': 'mecA'}])
    assert mec == ['mecA']
    assert score == 1
    assert special == set()


def test_amr_last_resort_resistance():
    mec, score, special = calculate_amr([{'gene': 'mecA'}, {'gene': 'vanA'}])
    assert mec == ['mecA']
    assert score == 2
    assert special == {'VANCOMYCIN'}


def test_cli_version():
    result = CliRunner().invoke(main, ['--version'])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_cli_rejects_invalid_database(tmp_path, monkeypatch):
    import coccivr.configure as configure

    monkeypatch.setattr(configure, 'exe', {tool: '/bin/true' for tool in configure.TOOLS})
    query = tmp_path / 'genome.fna'
    query.write_text('>c\nACGT\n')
    result = CliRunner().invoke(
        main, ['-q', str(query), '-p', 'x', '-o', str(tmp_path / 'out'),
               '--db', str(tmp_path / 'missing')])
    assert result.exit_code != 0
    assert 'Not a CocciVR database directory' in result.output
