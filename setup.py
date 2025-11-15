# Imports:
from setuptools import setup, find_packages

setup(
    name="gravmagpy",
    # Version scheme is last date updated in YYYY.MM.DDv format
    # where 'v' may increment [a...z] for multiple releases on the same day
    version="2025.11.15a",
    description="Gravity and Magnetic (GravMag) Modeling and Visualization Python Package",
    author="Dany Waller",
    author_email="dany.c.waller@gmail.com",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    python_requires=">=3.13, <4",
)