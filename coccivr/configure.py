# SPDX-License-Identifier: GPL-3.0-or-later
import os
import sys
from shutil import which


# When CocciVR lives in a conda environment, its companion tools (blast,
# diamond, amrfinder) are usually in the same environment. Make them
# discoverable even when the environment was not activated.
_env_bin = os.path.join(sys.prefix, 'bin')
if os.path.isdir(_env_bin) and _env_bin not in os.environ.get('PATH', '').split(os.pathsep):
    os.environ['PATH'] = _env_bin + os.pathsep + os.environ.get('PATH', '')

TOOLS = ('makeblastdb', 'blastn', 'diamond', 'amrfinder')
exe = {tool: which(tool) for tool in TOOLS}

# File names expected inside the database directory. Callers refer to the keys,
# so a file can be renamed here without touching the rest of the code.
DB_FILES = dict(
    vfdb='VFDB_setB_pro.fas',
    scc_type_db='whole_cassette_SCCmec_database_EXTENDED_20171117.fasta',
)

scheme_info = dict(min_iden=0.65, min_frag=0.6, max_iden=0.95)

db_folder = None
db_list = {}


def _db_candidates():
    """Places where a database directory is looked for, in priority order."""
    env = os.environ.get('COCCIVR_DB')
    if env:
        yield os.path.expanduser(env)
    here = os.path.dirname(os.path.abspath(__file__))
    # conda package layout: $CONDA_PREFIX/share/coccivr/db
    yield os.path.join(sys.prefix, 'share', 'staphmix', 'db')
    # source checkout: <project>/db, or <repo>/db when the project is a
    # subdirectory of the repository
    yield os.path.join(os.path.dirname(here), 'db')
    yield os.path.join(os.path.dirname(os.path.dirname(here)), 'db')
    yield os.path.join(here, 'db')


def is_db_dir(path):
    return bool(path) and os.path.isfile(os.path.join(path, DB_FILES['vfdb']))


def missing_db_files(folder):
    return [name for name in DB_FILES.values()
            if not os.path.isfile(os.path.join(folder, name))]


def default_db_folder():
    for path in _db_candidates():
        if is_db_dir(path):
            return os.path.abspath(path)
    return None


def set_db_folder(path):
    global db_folder, db_list
    folder = os.path.abspath(os.path.expanduser(path))
    if not is_db_dir(folder):
        raise ValueError(
            "Not a CocciVR database directory: {0}\n"
            "It must contain {1}.".format(folder, DB_FILES['vfdb']))
    db_folder = folder
    db_list = {key: os.path.join(folder, name) for key, name in DB_FILES.items()}
    return db_folder


def db_path(key):
    if not db_list:
        folder = default_db_folder()
        if folder is None:
            raise RuntimeError(
                "CocciVR database not found. Point --db or $COCCIVR_DB at a "
                "directory containing {0}.".format(DB_FILES['vfdb']))
        set_db_folder(folder)
    return db_list[key]


def check_tools():
    missing = [tool for tool, path in exe.items() if not path]
    if missing:
        raise RuntimeError(
            "Missing required executable(s): {0}. Install them and make sure "
            "they are on PATH.".format(', '.join(missing)))
