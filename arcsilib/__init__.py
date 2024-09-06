"""
ARCSI - this file is needed to ensure it can be imported

See other source files for details
"""

import os
from distutils.version import LooseVersion

ARCSI_VERSION_MAJOR = 4
ARCSI_VERSION_MINOR = 0
ARCSI_VERSION_PATCH = 6

ARCSI_VERSION = f"{ARCSI_VERSION_MAJOR}.{ARCSI_VERSION_MINOR}.{ARCSI_VERSION_PATCH}"
ARCSI_VERSION_OBJ = LooseVersion(ARCSI_VERSION)
__version__ = ARCSI_VERSION

ARCSI_COPYRIGHT_YEAR = "2015"
ARCSI_COPYRIGHT_NAMES = "Pete Bunting, Dan Clewley"

ARCSI_SUPPORT_EMAIL = "rsgislib-support@googlegroups.com"

ARCSI_WEBSITE = "https://github.com/remotesensinginfo/arcsi"

ARCSI_SENSORS_LIST = ["lsmss", "lstm", "lsetm", "lsoli", "sen2"]
ARCSI_PRODUCTS_LIST = [
    "RAD",
    "SATURATE",
    "TOA",
    "CLOUDS",
    "CLEARSKY",
    "DDVAOT",
    "DOSAOT",
    "DOSAOTSGL",
    "SREF",
    "STDSREF",
    "DOS",
    "THERMAL",
    "TOPOSHADOW",
    "FOOTPRINT",
    "METADATA",
    "SHARP",
]
ARCSI_ARCHIVE_EXE_LIST = [
    ".tar.gz",
    ".tgz",
    ".TAR.GZ",
    ".TGZ",
    ".tar",
    ".TAR",
    ".zip",
    ".ZIP",
    ".tar.bz",
    ".TAR.BZ",
    ".tar.bz2",
    ".TAR.BZ2",
]
ARCSI_GDALFORMATS_LIST = ["KEA"]
ARCSI_CLOUD_METHODS_LIST = [
    "FMASK",
    "FMASK_DISP",
    "S2CLOUDLESS",
    "S2LESSFMSK",
    "S2LESSFMSKD",
    "LSMSK",
]

# Get install prefix
install_prefix = __file__[: __file__.find("lib")]

DEFAULT_ARCSI_AEROIMG_PATH = os.path.join(
    install_prefix, "share", "arcsi", "WorldAerosolParams.kea"
)
DEFAULT_ARCSI_ATMOSIMG_PATH = os.path.join(
    install_prefix, "share", "arcsi", "WorldAtmosphereParams.kea"
)

# Check files exit - set to none if they don't
if not os.path.isfile(DEFAULT_ARCSI_ATMOSIMG_PATH):
    DEFAULT_ARCSI_ATMOSIMG_PATH = None
if not os.path.isfile(DEFAULT_ARCSI_AEROIMG_PATH):
    DEFAULT_ARCSI_AEROIMG_PATH = None

if os.environ.get("RIOS_DFLT_DRIVER", None) == None:
    os.environ["RIOS_DFLT_DRIVER"] = "KEA"
