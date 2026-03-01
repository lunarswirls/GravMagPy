# Imports:
from setuptools import setup, find_packages

setup(
    name="gravmagpy",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    # Version scheme is last date updated in YYYY.MM.DDv format
    # where 'v' may increment [a...z] for multiple releases on the same day
    version="2026.03.01a",
    description="Gravity and Magnetic (GravMag) Modeling and Visualization Python Package",
    author="Dany Waller",
    author_email="dany.c.waller@gmail.com",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    python_requires=">=3.13, <4",
)