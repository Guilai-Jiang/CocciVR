# Database licensing

The CocciVR code is licensed under GPL-3.0-or-later (see `LICENSE` and
`COPYING`). The reference data are separate works with their own terms. Before
redistributing the database bundle, confirm the terms below.

| File | Origin | When publishing |
|---|---|---|
| `VFDB_setB_pro.fas` | VFDB set B (Liu et al., *Nucleic Acids Research* 2022). | VFDB is provided for academic use; redistribution may require permission. Contact the VFDB maintainers before bundling. |
| `whole_cassette_SCCmec_database_EXTENDED_20171117.fasta` | Extended SCCmec cassette database used by SCCmecFinder (Kaya et al., *mSphere* 2018). | No explicit licence was found in the upstream repository at the time of writing. Obtain permission from the authors before redistributing. |
| AMRFinderPlus database | Downloaded separately by `amrfinder -u`. | Not bundled; governed by NCBI terms. |

If redistribution is not permitted, publish the code only and let users supply
their own database directory with `--db`.

## Citation

- Liu B, Zheng D, Zhou S, Chen L, Yang J. VFDB 2022: a general classification
  scheme for bacterial virulence factors. *Nucleic Acids Res*.
  2022;50(D1):D912–D917.
- Kaya H, et al. SCCmecFinder, a web-based tool for typing of staphylococcal
  cassette chromosome mec in Staphylococcus aureus using whole-genome sequence
  data. *mSphere*. 2018;3(1):e00612-17.
- Jiang G, et al. Large-scale genomic analysis places Chinese CC398 as a
  persistent human-associated MSSA lineage apart from the dominant global
  LA-MRSA clade. *mSystems*. 2026;11(10). doi:10.1128/msystems.00621-26
