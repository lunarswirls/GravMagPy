"""evaluate fortran volume fields at arbitrary orbital samples"""

import ctypes
from functools import lru_cache

import numpy as np

from .geometry import source_array
from .utils import build_fortran


@lru_cache(maxsize=8)
def load_backend(library):
    """load the c-interoperable double-precision fortran kernel"""
    backend = ctypes.CDLL(str(library))
    array = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    backend.orbital_field.argtypes = [
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double,
        array, array, array, array, array, ctypes.POINTER(ctypes.c_int),
    ]
    backend.orbital_field.restype = None
    return backend


@lru_cache(maxsize=8)
def load_dipole_backend(library):
    """load the dipole entrypoint without changing older volume-library requirements"""
    backend = ctypes.CDLL(str(library))
    array = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    backend.dipole_field.argtypes = [ctypes.c_int, ctypes.c_int, array, array, array, array,
                                    ctypes.POINTER(ctypes.c_int)]
    backend.dipole_field.restype = None
    return backend


@lru_cache(maxsize=16)
def quadrature(order):
    """cache gauss-legendre volume integration nodes and weights"""
    return tuple(np.ascontiguousarray(values) for values in np.polynomial.legendre.leggauss(order))


def predict_field(xyz_km, sources, *, radius_km=1737.4, quadrature_order=6, library=None):
    """return (n, 3) Bx/By/Bz in nt at actual planet-fixed positions in km

    source magnetization is in global planet-fixed xyz, not local north/east/down
    observers must be outside the reference sphere; volume quadrature should be
    increased to check convergence, especially near the surface
    """
    rows = source_array(sources, radius_km)
    xyz = np.ascontiguousarray(xyz_km, dtype=np.float64)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError("xyz_km must be a finite nonempty (n, 3) array")
    if np.any(np.linalg.norm(xyz, axis=1) <= radius_km):
        raise ValueError("orbital observations must be above the reference sphere")
    if not isinstance(quadrature_order, (int, np.integer)) or not 1 <= quadrature_order <= 64:
        raise ValueError("quadrature_order must be an integer between 1 and 64")
    if max(len(xyz), len(rows)) > np.iinfo(np.int32).max:
        raise ValueError("array dimensions exceed the fortran integer range")
    nodes, weights = quadrature(quadrature_order)
    backend = load_backend(library or build_fortran("orbital"))
    field = np.zeros_like(xyz)
    status = ctypes.c_int()
    backend.orbital_field(len(xyz), len(rows), quadrature_order, radius_km, xyz, rows, nodes, weights, field, ctypes.byref(status))
    if status.value or not np.isfinite(field).all():
        raise RuntimeError(f"Fortran orbital field failed with status {status.value}")
    return field
