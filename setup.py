"""include portable fortran sources in wheels without compiling at install time"""

from pathlib import Path
import shutil

from setuptools import setup
from setuptools.command.build_py import build_py


class build_with_fortran(build_py):
    def run(self):
        super().run()
        source = Path(__file__).parent / "fortran"
        destination = Path(self.build_lib) / "gravmagpy" / "_fortran"
        # rebuild the generated package source subtree to avoid stale files
        if destination.is_dir():
            shutil.rmtree(destination)
        shutil.copytree(source, destination, dirs_exist_ok=True, ignore=shutil.ignore_patterns("*.mod", "*.o"))


setup(cmdclass={"build_py": build_with_fortran})
