#
#  Java Information Dynamics Toolkit (JIDT)
#  Copyright (C) 2026, Alexander Richter
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""Places infodynamics.jar inside the package as the wheel is built.

The jar is taken from whichever source is available: one already in place,
one built from this checkout with ant (as per the AntScripts wiki page), or
the released full distribution for this version.
"""

import re
import shutil
import subprocess
import urllib.request
import zipfile
from pathlib import Path

from hatchling.builders.hooks.plugin.interface import BuildHookInterface
from hatchling.metadata.plugin.interface import MetadataHookInterface

HERE = Path(__file__).parent
REPO = HERE.parent
JAR = HERE / "src" / "jidt" / "infodynamics.jar"
DIST_URL = (
    "https://github.com/jlizier/jidt/releases/download"
    "/v{version}/infodynamics-dist-{version}.zip"
)


def toolkit_version() -> str:
    """Returns the version of the toolkit, as ant is told it in build.xml.

    build.xml holds the one version for the whole project, and the readme,
    the version file and the clojure project are all generated from it; the
    Python module is given it the same way rather than repeating it. An
    sdist carries no checkout, so there the version is read back from the
    PKG-INFO written when the sdist was built.

    Outputs:
    - result - the version, e.g. "1.6.1"
    """
    build = REPO / "build.xml"
    if build.exists():
        found = re.search(r'name="version" value="([^"]+)"', build.read_text())
        return found.group(1)
    # There is no checkout to read when building from an sdist, so take
    # the version recorded in it at the time it was made:
    pkg_info = HERE / "PKG-INFO"
    if not pkg_info.exists():
        raise FileNotFoundError(
            "Cannot tell which version of the toolkit this is: build from "
            "the python directory of a JIDT checkout, or from an sdist."
        )
    found = re.search(r"^Version: (.+)$", pkg_info.read_text(), re.MULTILINE)
    return found.group(1)


class VersionMetadataHook(MetadataHookInterface):
    """Supplies the toolkit's version as the version of this package."""

    PLUGIN_NAME = "custom"

    def update(self, metadata):
        metadata["version"] = toolkit_version()


def ensure_jar() -> Path:
    """Returns the jar for the package, fetching or building it if need be.

    Outputs:
    - result - path to infodynamics.jar inside the package
    """
    version = toolkit_version()
    if JAR.exists():
        return JAR
    if (REPO / "build.xml").exists() and shutil.which("ant"):
        try:
            subprocess.run(["ant", "jar"], cwd=REPO, check=True)
        except subprocess.CalledProcessError:
            # No Java development kit here, so fall back to the release:
            pass
        else:
            shutil.copyfile(REPO / "infodynamics.jar", JAR)
            return JAR
    # The full distribution for this release, as attached to its tag:
    zip_path = HERE / f"infodynamics-dist-{version}.zip"
    with (
        urllib.request.urlopen(DIST_URL.format(version=version)) as response,
        zip_path.open("wb") as archive,
    ):
        shutil.copyfileobj(response, archive)
    with zipfile.ZipFile(zip_path) as z, JAR.open("wb") as out:
        out.write(z.read("infodynamics.jar"))
    zip_path.unlink()
    return JAR


class JarBuildHook(BuildHookInterface):
    PLUGIN_NAME = "custom"

    def initialize(self, version, build_data):
        ensure_jar()
        build_data["artifacts"].append("src/jidt/infodynamics.jar")
