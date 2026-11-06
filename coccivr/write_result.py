# SPDX-License-Identifier: GPL-3.0-or-later
"""Output writers for the CocciVR program."""

import os


def write_scc(outdir_abpath, prefix, scc_out):
    output_name = prefix + '.sccmec.txt'
    header = [
        "Query contig", "Match ID", "Identity (%)", "Alignment Length", "Mismatch",
        "Strand", "Start Position (query)", "End Position (query)",
        "Start Position (subject)", "End Position (subject)", "E-value",
        "Bit Score", "Query Length", "Subject Length"]
    outfile = os.path.join(outdir_abpath, output_name)
    with open(outfile, 'wt') as fout:
        fout.write('\t'.join(header) + '\n')
        for value in scc_out.values():
            fout.write('\t'.join(str(item) for item in value[3][0]) + '\n')


def write_vf(outdir_abpath, prefix, vfdb_out):
    output_name = prefix + '.vf.txt'
    header = vfdb_out[0].keys()
    outfile = os.path.join(outdir_abpath, output_name)
    with open(outfile, 'wt') as fout:
        fout.write('\t'.join(header) + '\n')
        for entry in vfdb_out:
            fout.write('\t'.join(str(entry[key]) for key in header) + '\n')


def write_amr(outdir_abpath, prefix, amr_out):
    output_name = prefix + '.amr.txt'
    header = amr_out[0].keys()
    outfile = os.path.join(outdir_abpath, output_name)
    with open(outfile, 'wt') as fout:
        fout.write('\t'.join(header) + '\n')
        for entry in amr_out:
            fout.write('\t'.join(str(entry[key]) for key in header) + '\n')


def write_summary(outdir_abpath, prefix, summary_list):
    output_name = prefix + '.summary.txt'
    header = ['Query genome', 'SCCmec type', 'Enterotoxin score',
              'Blood infection score', 'MRSA/MSSA', 'Resistance score',
              'Antibiotic resistance gene_class',
              'Antibiotic resistance gene num', 'Special antibiotic resistance']
    outfile = os.path.join(outdir_abpath, output_name)
    with open(outfile, 'wt') as fout:
        fout.write('\t'.join(header) + '\n')
        fout.write('\t'.join(str(value) for value in summary_list) + '\n')
