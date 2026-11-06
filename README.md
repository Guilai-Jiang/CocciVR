# CocciVR

[![build](https://github.com/Guilai-Jiang/CocciVR/actions/workflows/build-conda.yml/badge.svg)](https://github.com/Guilai-Jiang/CocciVR/actions/workflows/build-conda.yml)
[![license](https://img.shields.io/badge/license-GPL--3.0--or--later-blue.svg)](LICENSE)

Virulence, resistance and SCCmec profiling for Gram-positive cocci. One
command runs the three modules and produces a single summary row.

## Requirements

Python ≥ 3.9 and, on `PATH`:

| Tool | Version |
|---|---|
| `makeblastdb`, `blastn` | ≥ 2.14 |
| `diamond` | ≥ 2.2 |
| `amrfinder` | ≥ 3.11 |

## Install

```bash
conda create -n coccivr -c conda-forge -c bioconda \
  python=3.10 blast diamond ncbi-amrfinderplus click
conda activate coccivr
pip install .
```

## Database

The database is distributed as a release asset:

```bash
curl -LO https://github.com/Guilai-Jiang/CocciVR/releases/download/db-v1.0/coccivr-db-1.0.tar.gz
mkdir -p ~/coccivr-db && tar xzf coccivr-db-1.0.tar.gz -C ~/coccivr-db
export COCCIVR_DB=~/coccivr-db
```

It is located in this order: `--db` → `$COCCIVR_DB` →
`$CONDA_PREFIX/share/coccivr/db` → `./db`.

| File | Content |
|---|---|
| `VFDB_setB_pro.fas` | VFDB set B proteins |
| `whole_cassette_SCCmec_database_EXTENDED_20171117.fasta` | SCCmec cassette database |

AMRFinderPlus downloads its own resistance database (`amrfinder -u`, or set
`AMRFINDER_DB`).

## Usage

```bash
coccivr -q genome.fna.gz -p sample -o results -n 8
```

| Option | Meaning |
|---|---|
| `-q, --query` | input genome: fasta/fastq, optionally gzipped |
| `-m, --sccmec` | run SCCmec typing |
| `-v, --vf` | run virulence gene comparison |
| `-a, --amr` | run resistance gene comparison |
| `-p, --prefix` | output file prefix |
| `-o, --outdir` | output directory |
| `-n, --n_thread` | threads (default 8) |
| `--db` | database directory |

With no module flag, all three modules run.

## Output

| File | Content |
|---|---|
| `*.summary.txt` | SCCmec type, virulence scores, MRSA/MSSA, resistance counts |
| `*.sccmec.txt` | SCCmec BLAST hits |
| `*.vf.txt` | virulence gene hits |
| `*.amr.txt` | resistance gene hits |

## Citation

> Jiang G, Guo C, Hua Y, Wu C, Li S, Wang Y, Zhao J, Ji L, Cheng Y, Zhou Z,
> Wu X, Li H. **Large-scale genomic analysis places Chinese CC398 as a
> persistent human-associated MSSA lineage apart from the dominant global
> LA-MRSA clade.** *mSystems*. 2026;11(10).
> [doi:10.1128/msystems.00621-26](https://doi.org/10.1128/msystems.00621-26)

## License

GPL-3.0-or-later. See [LICENSE](LICENSE) and [COPYING](COPYING). Database terms
are described in [docs/DATABASE-LICENSES.md](docs/DATABASE-LICENSES.md).

## Contact

Questions, bug reports and suggestions are very welcome. Please open an
[issue](https://github.com/Guilai-Jiang/CocciVR/issues) or write to
[gljiang@alu.suda.edu.cn](mailto:gljiang@alu.suda.edu.cn).
