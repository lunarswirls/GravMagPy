# GravMagPy
Gravity and magnetic equivalent source solver with minimization functions, global and local georeferenced map outputs, and field analysis and visualization tools.

<img width="2600" height="1250" alt="Golitsyn_context" src="https://github.com/user-attachments/assets/b84866f8-1d09-4200-a8ba-99e6abe7f0e6" />


## GravMag Sphere
Fortran code is provided in `gravmag_sphere` with [documentation](https://github.com/lunarswirls/GravMagPy/blob/73e48e9c7f4809ac19a8636c23b8d75ec3037ddf/gravmagpy/gravmag_sphere/README.md) and bash scripts to compile and run [examples](https://github.com/lunarswirls/GravMagPy/blob/73e48e9c7f4809ac19a8636c23b8d75ec3037ddf/gravmagpy/gravmag_sphere/examples/README_orientation_multibody_examples.md).

GravMag Sphere based on legacy `sphere` program developed by Dr. Dhahanjay Ravat, modernized and modularized during his Fall 2018 EES-395 "Gravity and Magnetics" course at the University of Kentucky.


## GravMagPy Installation
GravMagPy can be installed as a Python package called "gravmagpy" for easy access to functions and visualization tools. 
This was done so that one could organize files into separate directories and still be able to reference each file
by defining a path relative to the package root name.

To install as a local branch (highly recommended. NOTE the "-e" flag to install as a local branch! Again, highly 
recommended you do not ignore the "-e" flag!)

```
pip install -e <path/to/your/local/python/package/directory containing setup.py>
```
NOTE: This is not a path to setup.py! This is the path to the root directory containing setup.py!

As an example, if you are user tycho and you want to use your copy of the gravmagpy package that you put in 
`/homes/tycho/code/GravMagPy/gravmagpy`, you would do:
```
pip install -e /homes/tycho/code/GravMagPy/
```

## Fortran compilation
A bash script has been provided that compiles and runs ```gravmag_sphere_brtp```

Usage:
```
./run_gravmag_sphere.sh <R_sphere_km> <gauss_coeff_file> <input_in_file> [output_file]
```

An example input file has been provided that can be run as:
```
./run_gravmag_sphere.sh 1737.4 gauss2_32.data gravmag_sphere_test.in results.txt
```
and compared against the provided output file.

## Visuals
TBD


## Authors and acknowledgment
Primary Author: [Dany Waller](danywaller.github.io)
- Email: [dany.c.waller@gmail.com](mailto:dany.c.waller@gmail.com)
