"""
Routines for reading and writing PSI-style HDF5 and HDF4 data files.

This module provides a unified interface for interacting with PSI's HDF data
ecosystem.  It handles both HDF4 (``.hdf``) and HDF5 (``.h5``) file formats,
automatically dispatching to the appropriate backend based on the file
extension.

Key interfaces
--------------
Reading full datasets:
    :func:`read_hdf_data`, :func:`rdhdf_1d`, :func:`rdhdf_2d`, :func:`rdhdf_3d`

Writing full datasets:
    :func:`write_hdf_data`, :func:`wrhdf_1d`, :func:`wrhdf_2d`, :func:`wrhdf_3d`

Reading file metadata:
    :func:`read_hdf_meta`, :func:`read_rtp_meta`

Reading dataset subsets:
    :func:`get_scales_1d`, :func:`get_scales_2d`, :func:`get_scales_3d`,
    :func:`read_hdf_by_index`, :func:`read_hdf_by_value`, :func:`read_hdf_by_ivalue`

Interpolating data:
    :func:`np_interpolate_slice_from_hdf`, :func:`sp_interpolate_slice_from_hdf`,
    :func:`interpolate_positions_from_hdf`

Converting between formats:
    :func:`convert`, :func:`convert_psih4_to_psih5`

See Also
--------
:mod:`psi_data` :
    The ``psi-data-utils`` package, which provides helpers for fetching
    example HDF data files. See its `documentation
    <https://predsci.com/doc/psi-data-utils/index.html>`_ for details.

Examples
--------
Read a 3D PSI-style HDF5 file:

>>> from psi_io import read_hdf_data
>>> from psi_data import fetch_mas_data
>>> filepath = fetch_mas_data().cor_br
>>> data, r, t, p = read_hdf_data(filepath)
>>> data.shape
(299, 142, 255)
"""

from __future__ import annotations

__all__ = [
    "read_hdf_meta",
    "read_rtp_meta",

    "get_scales_1d",
    "get_scales_2d",
    "get_scales_3d",

    "read_hdf_by_index",
    "read_hdf_by_value",
    "read_hdf_by_ivalue",

    "np_interpolate_slice_from_hdf",
    "sp_interpolate_slice_from_hdf",
    "interpolate_positions_from_hdf",

    "instantiate_linear_interpolator",
    "interpolate_point_from_1d_slice",
    "interpolate_point_from_2d_slice",

    "read_hdf_data",
    "rdhdf_1d",
    "rdhdf_2d",
    "rdhdf_3d",

    "write_hdf_data",
    "write_hdf_meta",
    "wrhdf_1d",
    "wrhdf_2d",
    "wrhdf_3d",

    "convert",
    "convert_psih4_to_psih5"
]

import math
from collections import namedtuple
from pathlib import Path
from types import MappingProxyType
from typing import Optional, Literal, Tuple, Sequence, List, Dict, Union, Callable, Any, Mapping

import numpy as np
import h5py as h5

# -----------------------------------------------------------------------------
# Optional Imports and Import Checking
# -----------------------------------------------------------------------------
# These packages are needed by several functions and must be imported in the
# module namespace.
try:
    import pyhdf.SD as h4
    H4_AVAILABLE = True
    DTYPE_TO_SDC = MappingProxyType({
        "i": {
            1: h4.SDC.INT8,
            2: h4.SDC.INT16,
            4: h4.SDC.INT32,
        },
        "u": {
            1: h4.SDC.UINT8,
            2: h4.SDC.UINT16,
            4: h4.SDC.UINT32,
        },
        "f": {
            4: h4.SDC.FLOAT32,
            8: h4.SDC.FLOAT64,
        },
        "b": {
            1: h4.SDC.UINT8
        },
        "U": h4.SDC.CHAR,
        "S": h4.SDC.CHAR
    })
    """
        Helper dictionary mapping :class:`~numpy.dtype` kinds to HDF4 SDC types.

        The keys are :attr:`~numpy.dtype.kind`, and the values are either a direct
        mapping (for byte-string or unicode-string types) or a nested mapping of
        :attr:`~numpy.dtype.itemsize` to SDC type (for numeric types).
    """
except ImportError:
    H4_AVAILABLE = False
    DTYPE_TO_SDC = {}

try:
    from scipy.interpolate import RegularGridInterpolator
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


# Functions to stop execution if a package doesn't exist.
def _except_no_pyhdf():
    """Raise :exc:`ImportError` if the ``pyhdf`` package is not available.

    Raises
    ------
    ImportError
        If ``pyhdf`` was not importable at module load time.

    Examples
    --------
    >>> from psi_io.psi_io import _except_no_pyhdf, H4_AVAILABLE
    >>> H4_AVAILABLE  # doctest: +SKIP
    True
    >>> _except_no_pyhdf() is None  # no-op when pyhdf is installed
    True
    """
    if not H4_AVAILABLE:
        raise ImportError('The pyhdf package is required to read/write HDF4 .hdf files!')
    return


def _except_no_scipy():
    """Raise :exc:`ImportError` if the ``scipy`` package is not available.

    Raises
    ------
    ImportError
        If ``scipy`` was not importable at module load time.

    Examples
    --------
    >>> from psi_io.psi_io import _except_no_scipy, SCIPY_AVAILABLE
    >>> SCIPY_AVAILABLE  # doctest: +SKIP
    True
    >>> _except_no_scipy() is None  # no-op when scipy is installed
    True
    """
    if not SCIPY_AVAILABLE:
        raise ImportError('The scipy package is required for the interpolation routines!')
    return


SDC_TYPE_CONVERSIONS = MappingProxyType({
    3: np.dtype("ubyte"),
    4: np.dtype("byte"),
    5: np.dtype("float32"),
    6: np.dtype("float64"),
    20: np.dtype("int8"),
    21: np.dtype("uint8"),
    22: np.dtype("int16"),
    23: np.dtype("uint16"),
    24: np.dtype("int32"),
    25: np.dtype("uint32")
})
"""Helper dictionary for mapping HDF4 types to numpy dtypes"""


PSI_DATA_ID = MappingProxyType({
    'h4': 'Data-Set-2',
    'h5': 'Data'
})
"""Mapping of PSI standard dataset names for HDF4 and HDF5 files"""


PSI_SCALE_ID = MappingProxyType({
    'h4': ('fakeDim0', 'fakeDim1', 'fakeDim2'),
    'h5': ('dim1', 'dim2', 'dim3')
})
"""Mapping of PSI standard scale names for HDF4 and HDF5 files"""


HDFEXT = {'.hdf', '.h5'}
"""Set of possible HDF file extensions"""


HdfExtType = Literal['.hdf', '.h5']
"""Type alias for possible HDF file extensions"""


HdfScaleMeta = namedtuple('HdfScaleMeta', ['name', 'type', 'shape', 'attr', 'imin', 'imax'])
"""
    Named tuple storing metadata for a single HDF scale (coordinate) dimension.

    Parameters
    ----------
    name : str
        The name of the scale dataset.
    type : str
        The data type of the scale.
    shape : tuple[int, ...]
        The shape of the scale array.
    attr : dict[str, Any]
        A dictionary of attributes associated with the dataset.
    imin : float
        The minimum value of the scale.
        This assumes the scale is monotonically increasing.
    imax : float
        The maximum value of the scale.
        This assumes the scale is monotonically increasing.
"""


HdfDataMeta = namedtuple('HdfDataMeta', ['name', 'type', 'shape', 'attr', 'scales'])
"""
    Named tuple for HDF dataset metadata

    Parameters
    ----------
    name : str
        The name of the dataset.
    type : str
        The data type of the dataset.
    shape : tuple[int, ...]
        The shape of the dataset.
    attr : dict[str, Any]
        A dictionary of attributes associated with the dataset.
    scales : list[HdfScaleMeta]
        A list of scale metadata objects corresponding to each dimension of the dataset.
        If the dataset has no scales, this list will be empty.
"""

PathLike = Union[Path, str]
"""Type alias for file paths, accepting either :class:`pathlib.Path` or str"""


def _dtype_to_sdc(dtype: np.dtype):
    """Convert a numpy dtype to the corresponding HDF4 SDC type.

    Parameters
    ----------
    dtype : np.dtype
        The numpy dtype to convert.

    Returns
    -------
    out : int
        The corresponding HDF4 SDC type constant for the given numpy dtype.

    Raises
    ------
    ImportError
        If the ``pyhdf`` package is not available.
    KeyError
        If the dtype kind or itemsize is not supported by HDF4.  See
        :func:`write_hdf_data` for the full dtype support table.

    Examples
    --------
    >>> import numpy as np
    >>> from psi_io.psi_io import _dtype_to_sdc
    >>> import pyhdf.SD as h4
    >>> _dtype_to_sdc(np.dtype("float32")) == h4.SDC.FLOAT32
    True
    """
    _except_no_pyhdf()
    if dtype.kind in {"U", "S"}:
        return DTYPE_TO_SDC[dtype.kind]

    try:
        return DTYPE_TO_SDC[dtype.kind][dtype.itemsize]
    except KeyError as e:
        if dtype.kind not in DTYPE_TO_SDC:
            msg = (f"Unsupported dtype kind '{dtype.kind}' for HDF4. "
                   f"Supported kinds are: {set(DTYPE_TO_SDC.keys())}")
            raise KeyError(msg) from e
        elif dtype.itemsize not in DTYPE_TO_SDC[dtype.kind]:
            msg = (f"Unsupported itemsize '{dtype.itemsize}' for dtype kind '{dtype.kind}' in HDF4. "
                   f"Supported itemsizes are: {set(DTYPE_TO_SDC[dtype.kind].keys())}")
            raise KeyError(msg) from e
        raise e


def _dispatch_by_ext(ifile: PathLike,
                     hdf4_func: Callable,
                     hdf5_func: Callable,
                     *args: Any, **kwargs: Any
                     ):
    """
    Dispatch function to call HDF4 or HDF5 specific functions based on file extension.

    Parameters
    ----------
    ifile : PathLike
        The path to the HDF file.
    hdf4_func : Callable
        The function to call for HDF4 files.
    hdf5_func : Callable
        The function to call for HDF5 files.
    *args : Any
        Positional arguments to pass to the selected function.
    **kwargs : Any
        Keyword arguments to pass to the selected function.

    Returns
    -------
    out : Any
        The return value of the dispatched function.

    Raises
    ------
    ValueError
        If the file does not have a `.hdf` or `.h5` extension.
    ImportError
        If the file is HDF4 and the `pyhdf` package is not available.

    Examples
    --------
    >>> from psi_io.psi_io import _dispatch_by_ext, _read_h5_data, _read_h4_data
    >>> from psi_data import fetch_mas_data
    >>> filepath = fetch_mas_data().cor_br
    >>> data, *_ = _dispatch_by_ext(filepath, _read_h4_data, _read_h5_data)
    >>> data.shape
    (299, 142, 255)
    """
    ipath = Path(ifile)
    if ipath.suffix == '.h5':
        return hdf5_func(ifile, *args, **kwargs)
    if ipath.suffix == '.hdf':
        _except_no_pyhdf()
        return hdf4_func(ifile, *args, **kwargs)
    raise ValueError("File must be HDF4 (.hdf) or HDF5 (.h5)")


# -----------------------------------------------------------------------------
# "Classic" HDF reading and writing routines adapted from psihdf.py or psi_io.py.
# -----------------------------------------------------------------------------


def rdhdf_1d(hdf_filename: PathLike
             ) -> Tuple[np.ndarray, np.ndarray]:
    """Read a 1D PSI-style HDF5 or HDF4 file.

    Parameters
    ----------
    hdf_filename : PathLike
        The path to the 1D HDF5 (.h5) or HDF4 (.hdf) file to read.

    Returns
    -------
    x : np.ndarray
        1D coordinate scale array.
    f : np.ndarray
        1D data array.

    See Also
    --------
    read_hdf_data : Generic HDF data reading routine.

    Examples
    --------
    >>> from psi_io import rdhdf_1d
    >>> from psi_data import fetch_example_radial_scale
    >>> filepath = fetch_example_radial_scale()
    >>> x, f = rdhdf_1d(filepath)
    >>> x.shape, f.shape
    ((121,), (121,))
    """
    return _rdhdf_nd(hdf_filename, dimensionality=1)


def rdhdf_2d(hdf_filename: PathLike
             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a 2D PSI-style HDF5 or HDF4 file.

    The data in the HDF file is assumed to be ordered X, Y in Fortran order.
    Each dimension is assumed to have a 1D scale associated with it that
    describes the rectilinear grid coordinates.

    Parameters
    ----------
    hdf_filename : PathLike
        The path to the 2D HDF5 (.h5) or HDF4 (.hdf) file to read.

    Returns
    -------
    x : np.ndarray
        1D coordinate scale array for the X dimension.
    y : np.ndarray
        1D coordinate scale array for the Y dimension.
    f : np.ndarray
        2D data array, C-ordered with shape ``(ny, nx)`` in Python.

    Notes
    -----
    Because Fortran uses column-major ordering and Python uses row-major
    ordering, the data array is returned with reversed axis order relative to
    how it is stored on disk.  The first returned scale always corresponds to
    the *slowest* varying dimension in the data array (i.e. the last axis).

    See Also
    --------
    read_hdf_data : Generic HDF data reading routine.

    Examples
    --------
    >>> from psi_io import rdhdf_2d
    >>> from psi_data import fetch_mas_quantities
    >>> filepath = fetch_mas_quantities().ch_pm
    >>> x, y, f = rdhdf_2d(filepath)
    >>> f.shape == (y.shape[0], x.shape[0])
    True
    """
    return _rdhdf_nd(hdf_filename, dimensionality=2)


def rdhdf_3d(hdf_filename: PathLike
             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read a 3D PSI-style HDF5 or HDF4 file.

    The data in the HDF file is assumed to be ordered X, Y, Z in Fortran order.
    Each dimension is assumed to have a 1D scale associated with it that
    describes the rectilinear grid coordinates.

    Parameters
    ----------
    hdf_filename : PathLike
        The path to the 3D HDF5 (.h5) or HDF4 (.hdf) file to read.

    Returns
    -------
    x : np.ndarray
        1D coordinate scale array for the X dimension.
    y : np.ndarray
        1D coordinate scale array for the Y dimension.
    z : np.ndarray
        1D coordinate scale array for the Z dimension.
    f : np.ndarray
        3D data array, C-ordered with shape ``(nz, ny, nx)`` in Python.

    Notes
    -----
    Because Fortran uses column-major ordering and Python uses row-major
    ordering, the data array is returned with reversed axis order relative to
    how it is stored on disk.  For PSI spherical datasets the returned shapes
    follow the convention ``f.shape == (n_phi, n_theta, n_r)`` with scales
    returned in physical order ``(r, theta, phi)``.

    See Also
    --------
    read_hdf_data : Generic HDF data reading routine.

    Examples
    --------
    >>> from psi_io import rdhdf_3d
    >>> from psi_data import fetch_mas_data
    >>> filepath = fetch_mas_data().cor_br
    >>> r, t, p, f = rdhdf_3d(filepath)
    >>> f.shape == (p.shape[0], t.shape[0], r.shape[0])
    True
    """
    return _rdhdf_nd(hdf_filename, dimensionality=3)


def wrhdf_1d(hdf_filename: PathLike,
             x: np.ndarray,
             f: np.ndarray,
             **kwargs) -> None:
    """Write a 1D PSI-style HDF5 or HDF4 file.

    Parameters
    ----------
    hdf_filename : PathLike
        The path to the 1D HDF5 (.h5) or HDF4 (.hdf) file to write.
    x : np.ndarray
        1D array of scales.
    f : np.ndarray
        1D array of data.
    **kwargs
        Additional keyword arguments passed through to the underlying
        :func:`~psi_io.psi_io.write_hdf_data` routine, specifically:

        - ``dataset_id`` : str | None
            The identifier of the dataset to write.
            If None, a default dataset is used ('Data-Set-2' for HDF4 and 'Data' for HDF5).
        - ``sync_dtype``: bool
            If True, the data type of the scales will be matched to that of the data array.

        Omitting these will yield the same behavior as the legacy routines, *i.e.* writing to
        the default PSI dataset IDs for HDF4/HDF5 files and synchronizing datatypes between
        the dataset and scales.

    Raises
    ------
    ValueError
        If the file does not have a ``.hdf`` or ``.h5`` extension.
    KeyError
        If, for HDF4 files, the data or scale dtype is not supported by
        :py:mod:`pyhdf`.  See :func:`write_hdf_data` for the full dtype
        support table.

    Notes
    -----
    This routine is provided for backward compatibility with existing PSI codes that
    expect this API signature.  The underlying implementation dispatches to
    :func:`~psi_io.psi_io.write_hdf_data`, but the argument order and default
    behavior (sync dtype, default dataset IDs) are preserved.

    .. warning::
       When called with its default arguments this routine **writes to the default PSI
       dataset IDs for HDF4/HDF5 files and synchronizes dtypes between the dataset and
       scales.**  This behavior is required for interoperability with certain Fortran-based
       PSI tools.

    See Also
    --------
    write_hdf_data : Generic HDF data writing routine.

    Examples
    --------
    >>> import tempfile, numpy as np
    >>> from pathlib import Path
    >>> from psi_io import wrhdf_1d, rdhdf_1d
    >>> x = np.linspace(0.0, 1.0, 50, dtype=np.float32)
    >>> f = np.sin(x)
    >>> with tempfile.TemporaryDirectory() as d:
    ...     wrhdf_1d(Path(d) / "out.h5", x, f)
    ...     x2, f2 = rdhdf_1d(Path(d) / "out.h5")
    ...     np.array_equal(f, f2)
    True
    """
    return _wrhdf_nd(hdf_filename, f, x, dimensionality=1, **kwargs)


def wrhdf_2d(hdf_filename: PathLike,
             x: np.ndarray,
             y: np.ndarray,
             f: np.ndarray,
             **kwargs) -> None:
    """Write a 2D PSI-style HDF5 or HDF4 file.

    The data in the HDF file will appear as X,Y in Fortran order.

    Each dimension requires a 1D "scale" associated with it that
    describes the rectilinear grid coordinates in each dimension.

    Parameters
    ----------
    hdf_filename : PathLike
        The path to the 2D HDF5 (.h5) or HDF4 (.hdf) file to write.
    x : np.ndarray
        1D array of scales in the X dimension.
    y : np.ndarray
        1D array of scales in the Y dimension.
    f : np.ndarray
        2D data array, C-ordered with shape ``(ny, nx)`` in Python.
    **kwargs
        Additional keyword arguments passed through to the underlying
        :func:`~psi_io.psi_io.write_hdf_data` routine, specifically:

        - ``dataset_id`` : str | None
            The identifier of the dataset to write.
            If None, a default dataset is used (``'Data-Set-2'`` for HDF4 and ``'Data'`` for HDF5).
        - ``sync_dtype``: bool
            If True, the data type of the scales will be matched to that of the data array.

        Omitting these will yield the same behavior as the legacy routines, *i.e.* writing to
        the default PSI dataset IDs for HDF4/HDF5 files and synchronizing datatypes between
        the dataset and scales.

    Raises
    ------
    ValueError
        If the file does not have a ``.hdf`` or ``.h5`` extension.
    KeyError
        If, for HDF4 files, the data or scale dtype is not supported by
        :py:mod:`pyhdf`.  See :func:`write_hdf_data` for the full dtype
        support table.

    Notes
    -----
    This routine is provided for backward compatibility with existing PSI codes that
    expect this API signature.  The underlying implementation dispatches to
    :func:`~psi_io.psi_io.write_hdf_data`, but the argument order and default
    behavior (sync dtype, default dataset IDs) are preserved.

    .. warning::
       When called with its default arguments this routine **writes to the default PSI
       dataset IDs for HDF4/HDF5 files and synchronizes dtypes between the dataset and
       scales.**  This behavior is required for interoperability with certain Fortran-based
       PSI tools.

    See Also
    --------
    write_hdf_data : Generic HDF data writing routine.

    Examples
    --------
    >>> import tempfile, numpy as np
    >>> from pathlib import Path
    >>> from psi_io import wrhdf_2d, rdhdf_2d
    >>> x = np.linspace(0.0, 2*np.pi, 64, dtype=np.float32)
    >>> y = np.linspace(0.0, np.pi, 32, dtype=np.float32)
    >>> f = np.outer(np.sin(x), np.cos(y)).astype(np.float32)
    >>> with tempfile.TemporaryDirectory() as d:
    ...     wrhdf_2d(Path(d) / "out.h5", x, y, f)
    ...     x2, y2, f2 = rdhdf_2d(Path(d) / "out.h5")
    ...     np.array_equal(f, f2)
    True
    """
    return _wrhdf_nd(hdf_filename, f, x, y, dimensionality=2, **kwargs)


def wrhdf_3d(hdf_filename: PathLike,
             x: np.ndarray,
             y: np.ndarray,
             z: np.ndarray,
             f: np.ndarray,
             **kwargs) -> None:
    """Write a 3D PSI-style HDF5 or HDF4 file.

    The data in the HDF file will appear as X,Y,Z in Fortran order.

    Each dimension requires a 1D "scale" associated with it that
    describes the rectilinear grid coordinates in each dimension.

    Parameters
    ----------
    hdf_filename : PathLike
        The path to the 3D HDF5 (.h5) or HDF4 (.hdf) file to write.
    x : np.ndarray
        1D array of scales in the X dimension.
    y : np.ndarray
        1D array of scales in the Y dimension.
    z : np.ndarray
        1D array of scales in the Z dimension.
    f : np.ndarray
        3D data array, C-ordered with shape ``(nz, ny, nx)`` in Python.
    **kwargs
        Additional keyword arguments passed through to the underlying
        :func:`~psi_io.psi_io.write_hdf_data` routine, specifically:

        - ``dataset_id`` : str | None
            The identifier of the dataset to write.
            If None, a default dataset is used (``'Data-Set-2'`` for HDF4 and ``'Data'`` for HDF5).
        - ``sync_dtype``: bool
            If True, the data type of the scales will be matched to that of the data array.

        Omitting these will yield the same behavior as the legacy routines, *i.e.* writing to
        the default PSI dataset IDs for HDF4/HDF5 files and synchronizing datatypes between
        the dataset and scales.

    Raises
    ------
    ValueError
        If the file does not have a ``.hdf`` or ``.h5`` extension.
    KeyError
        If, for HDF4 files, the data or scale dtype is not supported by
        :py:mod:`pyhdf`.  See :func:`write_hdf_data` for the full dtype
        support table.

    Notes
    -----
    This routine is provided for backward compatibility with existing PSI codes that
    expect this API signature.  The underlying implementation dispatches to
    :func:`~psi_io.psi_io.write_hdf_data`, but the argument order and default
    behavior (sync dtype, default dataset IDs) are preserved.

    .. warning::
       When called with its default arguments this routine **writes to the default PSI
       dataset IDs for HDF4/HDF5 files and synchronizes dtypes between the dataset and
       scales.**  This behavior is required for interoperability with certain Fortran-based
       PSI tools.

    See Also
    --------
    write_hdf_data : Generic HDF data writing routine.

    Examples
    --------
    >>> import tempfile, numpy as np
    >>> from pathlib import Path
    >>> from psi_io import wrhdf_3d, rdhdf_3d
    >>> r = np.linspace(1.0, 5.0, 10, dtype=np.float32)
    >>> t = np.linspace(0.0, np.pi, 20, dtype=np.float32)
    >>> p = np.linspace(0.0, 2*np.pi, 30, dtype=np.float32)
    >>> f = np.ones((30, 20, 10), dtype=np.float32)
    >>> with tempfile.TemporaryDirectory() as d:
    ...     wrhdf_3d(Path(d) / "out.h5", r, t, p, f)
    ...     r2, t2, p2, f2 = rdhdf_3d(Path(d) / "out.h5")
    ...     np.array_equal(f, f2)
    True
    """
    return _wrhdf_nd(hdf_filename, f, x, y, z, dimensionality=3, **kwargs)


def get_scales_1d(filename: PathLike
                  ) -> np.ndarray:
    """Return the coordinate scale of a 1D PSI-style HDF5 or HDF4 dataset.

    Does not load the data array, so this is efficient for large files.

    Parameters
    ----------
    filename : PathLike
        The path to the 1D HDF5 (.h5) or HDF4 (.hdf) file to read.

    Returns
    -------
    x : np.ndarray
        1D coordinate scale array.

    Examples
    --------
    >>> from psi_io import get_scales_1d
    >>> from psi_data import fetch_example_radial_scale
    >>> filepath = fetch_example_radial_scale()
    >>> x = get_scales_1d(filepath)
    >>> x.ndim
    1
    """
    return _dispatch_by_ext(filename, _get_scales_nd_h4, _get_scales_nd_h5,
                            dimensionality=1)


def get_scales_2d(filename: PathLike
                  ) -> Tuple[np.ndarray, np.ndarray]:
    """Return the coordinate scales of a 2D PSI-style HDF5 or HDF4 dataset.

    Does not load the data array, so this is efficient for large files.
    The data in the HDF file is assumed to be ordered X, Y in Fortran order.

    Parameters
    ----------
    filename : PathLike
        The path to the 2D HDF5 (.h5) or HDF4 (.hdf) file to read.

    Returns
    -------
    x : np.ndarray
        1D coordinate scale array for the X dimension.
    y : np.ndarray
        1D coordinate scale array for the Y dimension.

    Examples
    --------
    >>> from psi_io import get_scales_2d
    >>> from psi_data import fetch_mas_quantities
    >>> filepath = fetch_mas_quantities().ch_pm
    >>> x, y = get_scales_2d(filepath)
    >>> x.ndim, y.ndim
    (1, 1)
    """
    return _dispatch_by_ext(filename, _get_scales_nd_h4, _get_scales_nd_h5,
                            dimensionality=2)


def get_scales_3d(filename: PathLike
                  ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the coordinate scales of a 3D PSI-style HDF5 or HDF4 dataset.

    Does not load the data array, so this is efficient for large files.
    The data in the HDF file is assumed to be ordered X, Y, Z in Fortran order.

    Parameters
    ----------
    filename : PathLike
        The path to the 3D HDF5 (.h5) or HDF4 (.hdf) file to read.

    Returns
    -------
    x : np.ndarray
        1D coordinate scale array for the X dimension.
    y : np.ndarray
        1D coordinate scale array for the Y dimension.
    z : np.ndarray
        1D coordinate scale array for the Z dimension.

    Examples
    --------
    >>> from psi_io import get_scales_3d
    >>> from psi_data import fetch_mas_data
    >>> filepath = fetch_mas_data().cor_br
    >>> r, t, p = get_scales_3d(filepath)
    >>> r.ndim, t.ndim, p.ndim
    (1, 1, 1)
    """
    return _dispatch_by_ext(filename, _get_scales_nd_h4, _get_scales_nd_h5,
                            dimensionality=3)


# -----------------------------------------------------------------------------
# "Updated" HDF reading and slicing routines for Hdf4 and Hdf5 datasets.
# -----------------------------------------------------------------------------


def read_hdf_meta(ifile: PathLike, /,
                  dataset_id: Optional[str] = None
                  ) -> List[HdfDataMeta]:
    """
    Read metadata from an HDF4 (.hdf) or HDF5 (.h5) file.

    This function provides a unified interface to read metadata from both HDF4 and HDF5 files.

    .. warning::
       Unlike elsewhere in this module, the scales and datasets are read **as is**, *i.e.* without
       reordering scales to match PSI's Fortran data ecosystem.

    .. warning::
       Unlike elsewhere in this module, when ``None`` is passed to ``dataset_id``, all (non-scale)
       datasets are returned (instead of the default psi datasets *e.g.* 'Data-Set-2' or 'Data').
       This will, effectively, return the standard PSI datasets when reading PSI-style HDF files.

    Parameters
    ----------
    ifile : PathLike
        The path to the HDF file to read.
    dataset_id : str | None, optional
        The identifier of the dataset for which to read metadata.
        If ``None``, metadata for **all** datasets is returned.  Default is ``None``.

    Returns
    -------
    out : list[HdfDataMeta]
        A list of :class:`HdfDataMeta` named tuples, one per dataset (or one
        for the requested *dataset_id*).  Each object exposes:

        ``name`` — dataset name (e.g. ``'Data'`` or ``'Data-Set-2'``).

        ``type`` — NumPy dtype string of the stored array.

        ``shape`` — array shape tuple.

        ``attr`` — :class:`dict` of attribute key-value pairs attached to the
        dataset.  Empty when no attributes are present.

        ``scales`` — list of :class:`HdfScaleMeta` objects, one per
        coordinate dimension.  Each :class:`HdfScaleMeta` similarly carries
        an ``attr`` dict for scale-level attributes.

    Raises
    ------
    ValueError
        If the file does not have a ``.hdf`` or ``.h5`` extension.

    Notes
    -----
    This function delegates to :func:`_read_h5_meta` for HDF5 files and
    :func:`_read_h4_meta` for HDF4 files based on the file extension.

    Although this function is designed to read metadata for dataset objects, it is possible
    to read metadata for coordinate variables (scales) by passing their names to
    ``dataset_id``, *e.g.* ``'dim1'``, ``'dim2'``, etc.  However, this is not the
    intended use case.

    Examples
    --------
    >>> from psi_io import read_hdf_meta
    >>> from psi_data import fetch_mas_data
    >>> filepath = fetch_mas_data().cor_br
    >>> meta = read_hdf_meta(filepath)
    >>> meta[0].name
    'Data'
    >>> meta[0].shape
    (299, 142, 255)
    >>> meta[0].attr                        # attributes attached to the dataset
    {}
    >>> len(meta[0].scales)                 # one HdfScaleMeta per dimension
    3
    """

    return _dispatch_by_ext(ifile, _read_h4_meta, _read_h5_meta,
                            dataset_id=dataset_id)


def read_rtp_meta(ifile: PathLike, /) -> Dict:
    """
    Read the scale metadata for PSI's 3D cubes.

    Parameters
    ----------
    ifile : PathLike
        The path to the HDF file to read.

    Returns
    -------
    out : dict
        A dictionary containing the RTP metadata.
        The value for each key ('r', 't', and 'p') is a tuple containing:

        1. The scale length
        2. The scale's value at the first index
        3. The scale's value at the last index

    Raises
    ------
    ValueError
        If the file does not have a `.hdf` or `.h5` extension.

    Notes
    -----
    This function delegates to :func:`_read_h5_rtp` for HDF5 files and
    :func:`_read_h4_rtp` for HDF4 files based on the file extension.

    Examples
    --------
    >>> from psi_io import read_rtp_meta
    >>> from psi_data import fetch_mas_data
    >>> filepath = fetch_mas_data().cor_br
    >>> meta = read_rtp_meta(filepath)
    >>> sorted(meta.keys())
    ['p', 'r', 't']
    >>> len(meta['r'])
    3
    """
    return _dispatch_by_ext(ifile, _read_h4_rtp, _read_h5_rtp)


def read_hdf_data(ifile: PathLike, /,
                  dataset_id: Optional[str] = None,
                  return_scales: bool = True,
                  ) -> Tuple[np.ndarray]:
    """
    Read data from an HDF4 (.hdf) or HDF5 (.h5) file.

    Parameters
    ----------
    ifile : PathLike
         The path to the HDF file to read.
    dataset_id : str | None, optional
        The identifier of the dataset to read.
        If ``None``, a default dataset is used (``'Data-Set-2'`` for HDF4 and
        ``'Data'`` for HDF5).  Default is ``None``.
    return_scales : bool, optional
        If ``True``, the coordinate scale arrays for each dimension are also
        returned.  Default is ``True``.

    Returns
    -------
    out : np.ndarray | tuple[np.ndarray, ...]
        The data array.  If ``return_scales`` is ``True``, returns a tuple
        ``(data, scale_0, scale_1, ...)`` with one scale per dimension.

    Raises
    ------
    ValueError
        If the file does not have a ``.hdf`` or ``.h5`` extension.

    See Also
    --------
    read_hdf_by_index : Read HDF datasets by index.
    read_hdf_by_value : Read HDF datasets by value ranges.
    read_hdf_by_ivalue : Read HDF datasets by subindex values.

    Notes
    -----
    This function delegates to :func:`_read_h5_data` for HDF5 files and
    :func:`_read_h4_data` for HDF4 files based on the file extension.

    Examples
    --------
    >>> from psi_io import read_hdf_data
    >>> from psi_data import fetch_mas_data
    >>> filepath = fetch_mas_data().cor_br
    >>> data, r, t, p = read_hdf_data(filepath)
    >>> data.shape
    (299, 142, 255)
    >>> r.shape, t.shape, p.shape
    ((255,), (142,), (299,))
    """
    return _dispatch_by_ext(ifile, _read_h4_data, _read_h5_data,
                            dataset_id=dataset_id, return_scales=return_scales)


def read_hdf_by_index(ifile: PathLike, /,
                      *xi: Union[int, Tuple[Union[int, None], Union[int, None]], None],
                      dataset_id: Optional[str] = None,
                      return_scales: bool = True,
                      ) -> Union[np.ndarray, Tuple[np.ndarray]]:
    r"""
    Read data from an HDF4 (.hdf) or HDF5 (.h5) file by index.

    .. attention::
       For each dimension, the *minimum* number of elements returned is 1 *e.g.*
       if 3 ints are passed (as positional `*xi` arguments) for a 3D dataset,
       the resulting subset will have a shape of (1, 1, 1,) with scales of length 1.

    Parameters
    ----------
    ifile : PathLike
       The path to the HDF file to read.
    *xi : int | tuple[int | None, int | None] | None
       Indices or ranges for each dimension of the `n`-dimensional dataset.
       Use None for a dimension to select all indices. If no arguments are passed,
       the entire dataset (and its scales) will be returned – see
       :func:`~psi_io.psi_io.read_hdf_data`.
    dataset_id : str | None, optional
       The identifier of the dataset to read.
       If ``None``, a default dataset is used (``'Data-Set-2'`` for HDF4 and
       ``'Data'`` for HDF5).  Default is ``None``.
    return_scales : bool, optional
       If ``True``, the coordinate scale arrays for each dimension are also
       returned.  Default is ``True``.

    Returns
    -------
    out : np.ndarray | tuple[np.ndarray, ...]
       The selected data array.  If ``return_scales`` is ``True``, returns a
       tuple ``(data, scale_0, scale_1, ...)`` with one scale per dimension.

    Raises
    ------
    ValueError
       If the file does not have a ``.hdf`` or ``.h5`` extension.

    See Also
    --------
    read_hdf_by_value : Read HDF datasets by value ranges.
    read_hdf_by_ivalue : Read HDF datasets by subindex values.
    read_hdf_data : Read entire HDF datasets.

    Notes
    -----
    This function delegates to :func:`_read_h5_by_index` for HDF5 files and
    :func:`_read_h4_by_index` for HDF4 files based on the file extension.

    This function assumes Fortran (column-major) ordering for compatibility with
    PSI's data ecosystem.  For an :math:`n`-dimensional array of shape
    :math:`(D_0, D_1, \ldots, D_{n-1})` the scales satisfy
    :math:`|x_i| = |D_{(n-1)-i}|`.  For example, a 3D dataset of shape
    :math:`(D_\phi, D_\theta, D_r)` has scales :math:`r, \theta, \phi`.

    Each ``*xi`` argument is forwarded to Python's built-in :class:`slice` to
    extract the desired subset without reading the entire dataset into memory.

    Examples
    --------
    Import a 3D HDF5 cube.

    >>> from psi_data import fetch_mas_data
    >>> from psi_io import read_hdf_by_index
    >>> filepath = fetch_mas_data().cor_br

    Extract a radial slice at the first radial index from a 3D cube:

    >>> f, r, t, p = read_hdf_by_index(filepath, 0, None, None)
    >>> f.shape, r.shape, t.shape, p.shape
    ((299, 142, 1), (1,), (142,), (299,))

    Extract a phi slice at the 90th index from a 3D cube:

    >>> f, r, t, p = read_hdf_by_index(filepath, None, None, 90)
    >>> f.shape, r.shape, t.shape, p.shape
    ((1, 142, 255), (255,), (142,), (1,))

    Extract up to the 20th radial index with phi indices 10 to 25:

    >>> f, r, t, p = read_hdf_by_index(filepath, (None, 20), None, (10, 25))
    >>> f.shape, r.shape, t.shape, p.shape
    ((15, 142, 20), (20,), (142,), (15,))
    """
    if not xi:
        return read_hdf_data(ifile, dataset_id=dataset_id, return_scales=return_scales)
    return _dispatch_by_ext(ifile, _read_h4_by_index, _read_h5_by_index,
                            *xi, dataset_id=dataset_id, return_scales=return_scales)


def read_hdf_by_value(ifile: PathLike, /,
                      *xi: Union[float, Tuple[float, float], None],
                      dataset_id: Optional[str] = None,
                      return_scales: bool = True,
                      ) -> Union[np.ndarray, Tuple[np.ndarray]]:
    """
    Read data from an HDF4 (.hdf) or HDF5 (.h5) file by value.

    .. note::
       For each dimension, the minimum number of elements returned is 2 *e.g.*
       if 3 floats are passed (as positional `*xi` arguments) for a 3D dataset,
       the resulting subset will have a shape of (2, 2, 2,) with scales of length 2.

    Parameters
    ----------
    ifile : PathLike
        The path to the HDF file to read.
    *xi : float | tuple[float, float] | None
        Values or value ranges corresponding to each dimension of the `n`-dimensional
        dataset specified by ``dataset_id``.  Pass ``None`` for a dimension to
        select all indices.  If no arguments are passed, the entire dataset (and
        its scales) will be returned.
    dataset_id : str | None, optional
        The identifier of the dataset to read.
        If ``None``, a default dataset is used (``'Data-Set-2'`` for HDF4 and
        ``'Data'`` for HDF5).  Default is ``None``.
    return_scales : bool, optional
        If ``True``, the coordinate scale arrays for each dimension are also
        returned.  Default is ``True``.

    Returns
    -------
    out : np.ndarray | tuple[np.ndarray, ...]
        The selected data array.  If ``return_scales`` is ``True``, returns a
        tuple ``(data, scale_0, scale_1, ...)`` with one scale per dimension.

    Raises
    ------
    ValueError
        If the file does not have a ``.hdf`` or ``.h5`` extension.

    See Also
    --------
    read_hdf_by_index : Read HDF datasets by index.
    read_hdf_by_ivalue : Read HDF datasets by subindex values.
    read_hdf_data : Read entire HDF datasets.
    sp_interpolate_slice_from_hdf : Interpolate slices using SciPy's
        :class:`~scipy.interpolate.RegularGridInterpolator`.
    np_interpolate_slice_from_hdf : Linear/bilinear/trilinear interpolation
        using vectorized NumPy routines.

    Notes
    -----
    This function delegates to :func:`_read_h5_by_value` for HDF5 files and
    :func:`_read_h4_by_value` for HDF4 files based on the file extension.

    This function assumes that the dataset is Fortran (or column-major) ordered *viz.* for
    compatibility with PSI's data ecosystem; as such, a given :math:`n`-dimensional array,
    of shape :math:`(D_0, D_1, ..., D_{n-1})`, has scales :math:`(x_0, x_1, ..., x_{n-1})`,
    such that :math:`| x_i | = | D_{(n-1)-i} |`. For example, a 3D dataset with shape
    :math:`(D_p, D_t, D_r)` has scales :math:`r, t, p` corresponding to the radial, theta,
    and phi dimensions respectively.

    This function extracts a subset of the given dataset/scales without reading the
    entire data into memory. For a given scale :math:`x_j`, if:

    - *i)* a single float is provided (:math:`a`), the function will return a 2-element
      subset of the scale (:math:`xʹ_j`) such that :math:`xʹ_j[0] <= a < xʹ_j[1]`.
    - *ii)* a (float, float) tuple is provided (:math:`a_0, a_1`), the function will return an
      *m*-element subset of the scale (:math:`xʹ_j`) where
      :math:`xʹ_j[0] <= a_0` and :math:`xʹ_j[m-1] > a_1`.
    - *iii)* a **None** value is provided, the function will return the entire scale :math:`x_j`

    The returned subset can then be passed to a linear interpolation routine to extract the
    "slice" at the desired fixed dimensions.

    Examples
    --------
    Import a 3D HDF5 cube.

    >>> from psi_data import fetch_mas_data
    >>> from psi_io import read_hdf_by_value
    >>> filepath = fetch_mas_data().cor_br

    Extract a radial slice at r=15 from a 3D cube:

    >>> f, r, t, p = read_hdf_by_value(filepath, 15, None, None)
    >>> f.shape, r.shape, t.shape, p.shape
    ((299, 142, 2), (2,), (142,), (299,))

    Extract a phi slice at p=1.57 from a 3D cube:

    >>> f, r, t, p = read_hdf_by_value(filepath, None, None, 1.57)
    >>> f.shape, r.shape, t.shape, p.shape
    ((2, 142, 255), (255,), (142,), (2,))

    Extract the values between 3.2 and 6.4 (in the radial dimension) and with
    phi equal to 4.5

    >>> f, r, t, p = read_hdf_by_value(filepath, (3.2, 6.4), None, 4.5)
    >>> f.shape, r.shape, t.shape, p.shape
    ((2, 142, 35), (35,), (142,), (2,))
    """
    if not xi:
        return read_hdf_data(ifile, dataset_id=dataset_id, return_scales=return_scales)
    return _dispatch_by_ext(ifile, _read_h4_by_value, _read_h5_by_value,
                            *xi, dataset_id=dataset_id, return_scales=return_scales)


def read_hdf_by_ivalue(ifile: PathLike, /,
                      *xi: Union[float, Tuple[float, float], None],
                      dataset_id: Optional[str] = None,
                      return_scales: bool = True,
                      ) -> Union[np.ndarray, Tuple[np.ndarray]]:
    r"""
    Read data from an HDF4 (.hdf) or HDF5 (.h5) file by subindex value.

    Unlike :func:`read_hdf_by_value` (which works in physical coordinate space),
    this function accepts fractional *index* positions and returns the minimal
    bracketing subset required for interpolation.

    .. note::
       For each dimension the minimum number of elements returned is 2.  For
       example, passing 3 scalar floats for a 3D dataset yields a subset of
       shape ``(2, 2, 2)`` with index arrays of length 2.

    Parameters
    ----------
    ifile : PathLike
        The path to the HDF file to read.
    *xi : float | tuple[float, float] | None
        Fractional index values or ranges for each dimension of the
        ``n``-dimensional dataset.  Pass ``None`` to select an entire dimension.
        If no arguments are passed, the full dataset is returned.
    dataset_id : str | None, optional
        The identifier of the dataset to read.  If ``None``, a default dataset
        is used (``'Data-Set-2'`` for HDF4 and ``'Data'`` for HDF5).
        Default is ``None``.
    return_scales : bool, optional
        If ``True``, 0-based index arrays (generated with :func:`~numpy.arange`)
        are returned for each dimension alongside the data.  These are always
        index-space arrays regardless of whether the dataset has physical
        coordinate scales.  Default is ``True``.

    Returns
    -------
    out : np.ndarray | tuple[np.ndarray, ...]
        The selected data array.  If ``return_scales`` is ``True``, returns a
        tuple ``(data, idx_0, idx_1, ...)`` with one index array per dimension.

    Raises
    ------
    ValueError
        If the file does not have a ``.hdf`` or ``.h5`` extension.

    See Also
    --------
    read_hdf_by_index : Read HDF datasets by integer index.
    read_hdf_by_value : Read HDF datasets by physical coordinate value.
    read_hdf_data : Read entire HDF datasets.
    sp_interpolate_slice_from_hdf : Interpolate slices using SciPy's
        :class:`~scipy.interpolate.RegularGridInterpolator`.
    np_interpolate_slice_from_hdf : Linear/bilinear/trilinear interpolation
        using vectorized NumPy routines.

    Notes
    -----
    This function delegates to :func:`_read_h5_by_ivalue` for HDF5 files and
    :func:`_read_h4_by_ivalue` for HDF4 files based on the file extension.

    For a given dimension with scale :math:`x_j`, the bracketing rule is:

    - *single float* :math:`a` → returns
      :math:`x_j[\lfloor a \rfloor], x_j[\lceil a \rceil]`.
    - *(float, float)* :math:`(a_0, a_1)` → returns
      :math:`x_j[\lfloor a_0 \rfloor], \ldots, x_j[\lceil a_1 \rceil]`.
    - ``None`` → returns the entire scale :math:`x_j`.

    Examples
    --------
    >>> from psi_data import fetch_mas_data
    >>> from psi_io import read_hdf_by_ivalue
    >>> filepath = fetch_mas_data().cor_br

    Extract a 2-element bracket around fractional radial index 2.7:

    >>> f, r_idx, t_idx, p_idx = read_hdf_by_ivalue(filepath, 2.7, None, None)
    >>> r_idx.shape
    (2,)
    """
    if not xi:
        return read_hdf_data(ifile, dataset_id=dataset_id, return_scales=return_scales)
    return _dispatch_by_ext(ifile, _read_h4_by_ivalue, _read_h5_by_ivalue,
                            *xi, dataset_id=dataset_id, return_scales=return_scales)


def write_hdf_data(ifile: PathLike, /,
                   data: np.ndarray,
                   *scales: Sequence[Union[np.ndarray, None]],
                   dataset_id: Optional[str] = None,
                   sync_dtype: bool = False,
                   strict: bool = True,
                   **kwargs
                   ) -> Path:
    """
    Write data to an HDF4 (.hdf) or HDF5 (.h5) file.

    Following PSI conventions, the data array is assumed to be Fortran-ordered,
    with the scales provided in the order corresponding to each dimension *e.g.* a
    3D dataset with shape (Dp, Dt, Dr) has scales r, t, p corresponding to the
    radial, theta, and phi dimensions respectively.

    Parameters
    ----------
    ifile : PathLike
        The path to the HDF file to write.
    data : np.ndarray
        The data array to write.
    *scales : Sequence[np.ndarray | None]
        The scales (coordinate arrays) for each dimension.
    dataset_id : str | None, optional
        The identifier of the dataset to write.  If ``None``, a default dataset
        is used (``'Data-Set-2'`` for HDF4 and ``'Data'`` for HDF5).
        Default is ``None``.
    sync_dtype : bool, optional
        If ``True``, the scale dtypes are cast to match the data array dtype.
        This mimics the behavior of PSI's legacy HDF writing routines and ensures
        compatibility with Fortran tools that require uniform precision between
        datasets and their scales.  Default is ``False``.
    strict : bool, optional
        If ``True``, raise an error if any dataset attribute cannot be written to
        the target format.  If ``False``, a warning is printed and the attribute
        is skipped.  Default is ``True``.
    **kwargs
        Key-value pairs of dataset attributes to attach to the dataset.

    Returns
    -------
    out : Path
        The path to the written HDF file.

    Raises
    ------
    ValueError
        If the file does not have a ``.hdf`` or ``.h5`` extension.
    KeyError
        If, for HDF4 files, the data or scale dtype is not supported by
        :py:mod:`pyhdf`.  See the dtype support table in the Notes section.

    Notes
    -----
    This function delegates to :func:`_write_h5_data` for HDF5 files and
    :func:`_write_h4_data` for HDF4 files based on the file extension.

    If no scales are provided the dataset is written without coordinate variables.
    The number of scales may be less than or equal to the number of dimensions;
    pass ``None`` for dimensions that should not have an attached scale.

    The table below summarizes dtype support across formats.  HDF4 support is
    determined by the :data:`DTYPE_TO_SDC` mapping; HDF5 support is provided
    by :mod:`h5py`.

    .. list-table:: NumPy dtype support by format
       :header-rows: 1
       :widths: 20 20 20

       * - dtype
         - HDF4 (``.hdf``)
         - HDF5 (``.h5``)
       * - ``float16``
         - No
         - Yes
       * - ``float32``
         - Yes
         - Yes
       * - ``float64``
         - Yes
         - Yes
       * - ``int8``
         - Yes
         - Yes
       * - ``int16``
         - Yes
         - Yes
       * - ``int32``
         - Yes
         - Yes
       * - ``int64``
         - No
         - Yes
       * - ``uint8``
         - Yes
         - Yes
       * - ``uint16``
         - Yes
         - Yes
       * - ``uint32``
         - Yes
         - Yes
       * - ``uint64``
         - No
         - Yes
       * - ``complex64``
         - No
         - Yes
       * - ``complex128``
         - No
         - Yes
       * - ``bool``
         - Yes (stored as ``uint8``)
         - Yes

    See Also
    --------
    wrhdf_1d : Write 1D HDF files.
    wrhdf_2d : Write 2D HDF files.
    wrhdf_3d : Write 3D HDF files.

    Examples
    --------
    >>> import tempfile, numpy as np
    >>> from pathlib import Path
    >>> from psi_io import write_hdf_data, read_hdf_data
    >>> r = np.linspace(1.0, 5.0, 10, dtype=np.float32)
    >>> t = np.linspace(0.0, np.pi, 20, dtype=np.float32)
    >>> p = np.linspace(0.0, 2*np.pi, 30, dtype=np.float32)
    >>> f = np.ones((30, 20, 10), dtype=np.float32)
    >>> with tempfile.TemporaryDirectory() as d:
    ...     _ = write_hdf_data(Path(d) / "out.h5", f, r, t, p)
    ...     data, r2, t2, p2 = read_hdf_data(Path(d) / "out.h5")
    ...     data.shape
    (30, 20, 10)
    """
    return _dispatch_by_ext(ifile, _write_h4_data, _write_h5_data, data,
                            *scales, dataset_id=dataset_id, sync_dtype=sync_dtype, strict=strict, **kwargs)


def write_hdf_meta(ifile: PathLike, /,
                   meta: Optional[Mapping[str, Mapping[str, Any]]] = None,
                   **kwargs) -> Path:
    """
    Add or update attributes on an existing HDF4 (.hdf) or HDF5 (.h5) file.

    Unlike :func:`write_hdf_data`, this function does **not** create or overwrite
    any dataset — it only attaches attribute key-value pairs to the file root
    (via ``**kwargs``) or to named datasets within the file (via ``meta``).

    Parameters
    ----------
    ifile : PathLike
        Path to an existing HDF file.  The file must already contain any
        datasets referenced by *meta*.
    meta : Mapping[str, Mapping[str, Any]] | None, optional
        Per-dataset attributes.  Each key is a dataset name (e.g. ``'Data'``
        for HDF5, ``'Data-Set-2'`` for HDF4); the corresponding value is a
        mapping of attribute name → attribute value to attach to that dataset.
        ``None`` (the default) applies no per-dataset attributes.
    **kwargs
        File-level (root) attributes.  For HDF5 files these are written to
        the root group; for HDF4 files they are written as SD-level global
        attributes.  They do **not** appear on any dataset and will not be
        returned by :func:`read_hdf_meta`.

    Returns
    -------
    out : Path
        The path to the modified HDF file.

    Raises
    ------
    ValueError
        If the file does not have a ``.hdf`` or ``.h5`` extension.

    Notes
    -----
    This function delegates to :func:`_write_h5_meta` for HDF5 files and
    :func:`_write_h4_meta` for HDF4 files based on the file extension.

    Per-dataset attributes written via *meta* are subsequently accessible
    through the :attr:`~HdfDataMeta.attr` field of the :class:`HdfDataMeta`
    named tuple returned by :func:`read_hdf_meta`.

    See Also
    --------
    write_hdf_data : Write dataset data (and optionally attributes) to a new file.
    read_hdf_meta  : Read dataset metadata and attributes from an HDF file.

    Examples
    --------
    >>> import tempfile, numpy as np
    >>> from pathlib import Path
    >>> from psi_io import write_hdf_data, write_hdf_meta, read_hdf_meta
    >>> f = np.ones((30, 20, 10), dtype=np.float32)
    >>> with tempfile.TemporaryDirectory() as d:
    ...     _ = write_hdf_data(Path(d) / "out.h5", f)
    ...     _ = write_hdf_meta(Path(d) / "out.h5", meta={'Data': {'long_name': 'test field', 'units': 'G'}})
    ...     meta = read_hdf_meta(Path(d) / "out.h5")
    ...     meta[0].attr['long_name']
    'test field'
    """
    return _dispatch_by_ext(ifile, _write_h4_meta, _write_h5_meta, meta, **kwargs)


def convert(ifile: PathLike,
            ofile: Optional[PathLike] = None,
            strict: bool = True) -> Path:
    """
    Convert an HDF file between HDF4 (.hdf) and HDF5 (.h5) formats.

    All datasets and their attributes are preserved in the output file.
    The output format is inferred from the input: `.hdf` → `.h5` and vice versa,
    unless an explicit output path is provided.

    Parameters
    ----------
    ifile : PathLike
        The path to the source HDF file.
    ofile : PathLike, optional
        The path to the output HDF file.
        If ``None``, the output file is written alongside the input file
        with its extension swapped (*e.g.* ``foo.hdf`` → ``foo.h5``).
    strict : bool, optional
        If ``True``, raise an error if any dataset attribute cannot be written
        to the output format.  If ``False``, a warning is printed and the
        attribute is skipped.  Default is ``True``.

    Returns
    -------
    out : Path
        The path to the written output file.

    Raises
    ------
    ValueError
        If the input file does not have a ``.hdf`` or ``.h5`` extension.
    KeyError
        If, for HDF4 output files, the data or an attribute value has a dtype
        not supported by :py:mod:`pyhdf` and ``strict`` is ``True``.  See
        :func:`write_hdf_data` for the full dtype support table.
    TypeError
        If, for HDF5 output files, a dataset attribute has a type that cannot be
        stored as an HDF5 attribute and ``strict`` is ``True``.

    See Also
    --------
    convert_psih4_to_psih5 : Specialized PSI-convention HDF4 → HDF5 converter.
    write_hdf_data : Generic HDF data writing routine.
    read_hdf_data : Generic HDF data reading routine.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> from psi_io import convert, read_hdf_meta
    >>> from psi_data import fetch_mas_data
    >>> h5_path = fetch_mas_data(hdf=5).cor_br
    >>> with tempfile.TemporaryDirectory() as d:
    ...     out = convert(h5_path, Path(d) / "br.hdf")
    ...     meta = read_hdf_meta(out)
    ...     meta[0].name
    'Data'
    """
    ifile = Path(ifile)
    if not ofile:
        ofile = ifile.with_suffix(".hdf") if ifile.suffix == ".h5" else ifile.with_suffix(".h5")
    else:
        ofile = Path(ofile)
    ofile.parent.mkdir(parents=True, exist_ok=True)

    meta_data = read_hdf_meta(ifile)
    for dataset in meta_data:
        data, *scales = read_hdf_data(ifile, dataset_id=dataset.name, return_scales=True)
        write_hdf_data(ofile, data, *scales, dataset_id=dataset.name, strict=strict, **dataset.attr)

    return ofile


def convert_psih4_to_psih5(ifile: PathLike,
                          ofile: Optional[PathLike] = None) -> Path:
    """
    Convert a PSI-convention HDF4 file to an HDF5 file.

    Unlike :func:`convert`, this function is specialized for PSI-style HDF4 files:
    it reads the primary dataset (``'Data-Set-2'``) and writes it under the PSI HDF5
    dataset name (``'Data'``), preserving the dataset's attributes and scales.

    Parameters
    ----------
    ifile : PathLike
        The path to the source HDF4 (.hdf) file.
    ofile : PathLike, optional
        The path to the output HDF5 (.h5) file.
        If ``None``, the output file is written alongside the input file with a
        ``.h5`` extension (*e.g.* ``foo.hdf`` → ``foo.h5``).

    Returns
    -------
    out : Path
        The path to the written HDF5 file.

    Raises
    ------
    ValueError
        If ``ifile`` does not have a ``.hdf`` extension.
    ValueError
        If ``ofile`` does not have a ``.h5`` extension.
    ImportError
        If the :py:mod:`pyhdf` package is not available.

    See Also
    --------
    convert : General-purpose HDF4 ↔ HDF5 converter.
    write_hdf_data : Generic HDF data writing routine.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> from psi_io import convert_psih4_to_psih5, read_hdf_meta
    >>> from psi_data import fetch_mas_data
    >>> hdf4_path = fetch_mas_data(hdf=4).cor_br
    >>> with tempfile.TemporaryDirectory() as d:
    ...     out = convert_psih4_to_psih5(hdf4_path, Path(d) / "br.h5")
    ...     meta = read_hdf_meta(out)
    ...     meta[0].name
    'Data'
    """
    ifile = Path(ifile)
    ofile = ifile.with_suffix(".h5") if not ofile else Path(ofile)
    if ifile.suffix != ".hdf":
        raise ValueError(f"Input file must have a .hdf extension; got {ifile.suffix}")
    if ofile.suffix != ".h5":
        raise ValueError(f"Output file must have a .h5 extension; got {ofile.suffix}")
    ofile.parent.mkdir(parents=True, exist_ok=True)


    data, *scales = read_hdf_data(ifile, dataset_id=PSI_DATA_ID["h4"], return_scales=True)
    meta_data, *_ = read_hdf_meta(ifile, dataset_id=PSI_DATA_ID["h4"])

    write_hdf_data(ofile, data, *scales,
                   dataset_id=PSI_DATA_ID["h5"], **meta_data.attr)
    return ofile

def instantiate_linear_interpolator(*args, **kwargs):
    r"""
    Instantiate a linear interpolator using the provided data and scales.

    Parameters
    ----------
    *args : sequence[array_like]
        The first argument is the data array.
        Subsequent arguments are the scales (coordinate arrays) for each dimension.
    **kwargs : dict
        Additional keyword arguments to pass to
        :class:`~scipy.interpolate.RegularGridInterpolator`.

    Returns
    -------
    out : RegularGridInterpolator
        An instance of RegularGridInterpolator initialized
        with the provided data and scales.

    Raises
    ------
    ImportError
        If the ``scipy`` package is not available.

    Notes
    -----
    This function transposes the data array and passes it along with the scales
    to :class:`~scipy.interpolate.RegularGridInterpolator`.  Given a PSI-style
    Fortran-ordered 3D dataset, the resulting interpolator can be queried using
    :math:`(r, \theta, \phi)` coordinates.

    Examples
    --------
    Import a 3D HDF5 cube.

    >>> from psi_data import fetch_mas_data
    >>> from psi_io import read_hdf_by_value
    >>> from numpy import pi
    >>> filepath = fetch_mas_data().cor_br

    Read the dataset by value (at 15 R_sun in the radial dimension).

    >>> data_and_scales = read_hdf_by_value(filepath, 15, None, None)
    >>> interpolator = instantiate_linear_interpolator(*data_and_scales)

    Interpolate at a specific position.

    >>> interpolator((15, pi/2, pi))
    array(-0.00117191)
    """
    _except_no_scipy()
    return RegularGridInterpolator(
        values=args[0].T,
        points=args[1:],
        **kwargs)


def sp_interpolate_slice_from_hdf(*xi, **kwargs):
    r"""
    Interpolate a slice from HDF data using SciPy's `RegularGridInterpolator`.

    .. note::
       Slicing routines result in a dimensional reduction. The dimensions
       that are fixed (i.e. provided as float values in `*xi`) are removed
       from the output slice, while the dimensions that are not fixed
       (*i.e.* provided as None in `*xi`) are retained.

    Parameters
    ----------
    *xi : sequence
        Positional arguments passed-through to :func:`read_hdf_by_value`.
    **kwargs : dict
        Keyword arguments passed-through to :func:`read_hdf_by_value`.
        **NOTE: Instantiating a linear interpolator requires the** ``return_scales``
        **keyword argument to be set to True; this function overrides
        any provided value for** ``return_scales`` **to ensure this behavior.**

    Returns
    -------
    data_slice : np.ndarray
        The interpolated data slice with fixed dimensions removed.
    scales : list[np.ndarray]
        Coordinate scale arrays for the retained (non-fixed) dimensions.

    Notes
    -----
    This function reads data from an HDF file, builds a
    :class:`~scipy.interpolate.RegularGridInterpolator`, and evaluates it at
    the requested fixed coordinates to produce the slice.

    .. note::
       The returned slice is Fortran-ordered, *e.g.* radial slices have shape
       :math:`(n_\phi, n_\theta)`, and :math:`\phi` slices have shape
       :math:`(n_r, n_\theta)`.

    .. note::
       :class:`~scipy.interpolate.RegularGridInterpolator` casts all input data
       to ``float64`` internally.  PSI HDF datasets stored as ``float32`` will
       therefore be upcast during interpolation.

    Examples
    --------
    >>> from psi_data import fetch_mas_data
    >>> from psi_io import sp_interpolate_slice_from_hdf
    >>> from numpy import pi
    >>> filepath = fetch_mas_data().cor_br

    Fetch a 2D slice at r=15 from 3D map

    >>> slice_, theta_scale, phi_scale = sp_interpolate_slice_from_hdf(filepath, 15, None, None)
    >>> slice_.shape, theta_scale.shape, phi_scale.shape
    ((299, 142), (142,), (299,))

    Fetch a single point from 3D map

    >>> point_value, *_ = sp_interpolate_slice_from_hdf(filepath, 1, pi/2, pi)
    >>> point_value
    np.float64(-0.20146621106799947)
    """
    filepath, *args = xi
    kwargs.pop('return_scales', None)
    result = read_hdf_by_value(filepath, *args, **kwargs)
    interpolator = instantiate_linear_interpolator(*result)
    grid = [yi[0] if yi[0] is not None else yi[1] for yi in zip(args, result[1:])]
    slice_ = interpolator(tuple(np.meshgrid(*grid, indexing='ij')))
    indices = [0 if yi is not None else slice(None, None) for yi in args]
    return slice_[tuple(indices)].T, *[yi[1] for yi in zip(args, result[1:]) if yi[0] is None]


def np_interpolate_slice_from_hdf(ifile: PathLike, /,
                       *xi: Union[float, Tuple[float, float], None],
                       dataset_id: Optional[str] = None,
                       by_index: bool = False,
                       ):
    """
    Interpolate a slice from HDF data using linear interpolation.

    .. note::
       Slicing routines result in a dimensional reduction. The dimensions
       that are fixed (i.e. provided as float values in `*xi`) are removed
       from the output slice, while the dimensions that are not fixed
       (*i.e.* provided as `None` in `*xi`) are retained.

    Parameters
    ----------
    ifile : PathLike
        The path to the HDF file to read.
    *xi : sequence
        Positional arguments passed-through to reader function.
    dataset_id : str | None, optional
        The identifier of the dataset to read.  If ``None``, a default dataset
        is used (``'Data-Set-2'`` for HDF4 and ``'Data'`` for HDF5).
        Default is ``None``.
    by_index : bool, optional
        If ``True``, use :func:`read_hdf_by_ivalue` to read data by fractional
        index values.  If ``False``, use :func:`read_hdf_by_value` to read by
        physical coordinate values.  Default is ``False``.

    Returns
    -------
    data_slice : np.ndarray
        The interpolated data slice with fixed dimensions removed.
    scales : list[np.ndarray]
        Coordinate scale arrays for the retained (non-fixed) dimensions.

    Raises
    ------
    ValueError
        If the number of dimensions to interpolate over is not supported.

    Notes
    -----
    This function supports linear, bilinear, and trilinear interpolation
    depending on the number of dimensions fixed in `xi`.

    Examples
    --------
    >>> from psi_data import fetch_mas_data
    >>> from psi_io import np_interpolate_slice_from_hdf
    >>> from numpy import pi
    >>> filepath = fetch_mas_data().cor_br

    Fetch a 2D slice at r=15 from 3D map

    >>> slice_, theta_scale, phi_scale = np_interpolate_slice_from_hdf(filepath, 15, None, None)
    >>> slice_.shape, theta_scale.shape, phi_scale.shape
    ((299, 142), (142,), (299,))

    Fetch a single point from 3D map

    >>> point_value, *_ = np_interpolate_slice_from_hdf(filepath, 1, pi/2, pi)
    >>> point_value
    np.float32(-0.2014672)

    """
    reader = read_hdf_by_value if not by_index else read_hdf_by_ivalue
    data, *scales = reader(ifile, *xi, dataset_id=dataset_id, return_scales=True)
    f_ = np.transpose(data)
    slice_type = sum([yi is not None for yi in xi])
    if slice_type == 1:
        return _np_linear_interpolation(xi, scales, f_).T, *[yi[1] for yi in zip(xi, scales) if yi[0] is None]
    elif slice_type == 2:
        return _np_bilinear_interpolation(xi, scales, f_).T, *[yi[1] for yi in zip(xi, scales) if yi[0] is None]
    elif slice_type == 3:
        return _np_trilinear_interpolation(xi, scales, f_).T, *[yi[1] for yi in zip(xi, scales) if yi[0] is None]
    else:
        raise ValueError("Not a valid number of dimensions for supported linear interpolation methods")


def interpolate_positions_from_hdf(ifile, *xi, **kwargs):
    r"""
    Interpolate at a list of scale positions using SciPy's
    :class:`~scipy.interpolate.RegularGridInterpolator`.

    Parameters
    ----------
    ifile : PathLike
        The path to the HDF file to read.
    *xi : sequence[np.ndarray]
       Coordinate values for each dimension of the ``n``-dimensional dataset.
       Each array must have the same length :math:`m`; the function assembles
       them into an :math:`m \times n` column stack for interpolation.
    **kwargs
        Keyword arguments forwarded to :func:`read_hdf_by_value`.

    Returns
    -------
    out : np.ndarray
        The interpolated values at the provided positions.

    Notes
    -----
    This function reads data from an HDF file, creates a linear interpolator,
    and interpolates at the provided scale values. For each dimension, the
    minimum and maximum values from the provided arrays are used to read
    the necessary subset of data from the HDF file *viz.* to avoid loading
    the entire dataset into memory.

    Examples
    --------
    Import a 3D HDF5 cube.

    >>> from psi_data import fetch_mas_data
    >>> from psi_io import interpolate_positions_from_hdf
    >>> import numpy as np
    >>> filepath = fetch_mas_data().cor_br

    Set up positions to interpolate.

    >>> r_vals = np.array([15, 20, 25])
    >>> theta_vals = np.array([np.pi/4, np.pi/2, 3*np.pi/4])
    >>> phi_vals = np.array([0, np.pi, 2*np.pi])

    Interpolate at the specified positions.

    >>> interpolate_positions_from_hdf(filepath, r_vals, theta_vals, phi_vals)
    array([-1.44383936e-03, -6.70081301e-04,  8.56632460e-05])
    """
    xi_ = [(np.nanmin(i), np.nanmax(i)) for i in xi]
    f, *scales = read_hdf_by_value(ifile, *xi_, **kwargs)
    interpolator = instantiate_linear_interpolator(f, *scales, bounds_error=False)
    return interpolator(np.stack(xi, axis=len(xi[0].shape)))


def interpolate_point_from_1d_slice(xi, scalex, values):
    """
    Interpolate a point from a 1D slice using linear interpolation.

    Parameters
    ----------
    xi : float
        The coordinate value at which to interpolate.
    scalex : np.ndarray
        The coordinate scale array for the dimension.
    values : np.ndarray
        The 1D data array to interpolate.

    Returns
    -------
    out : np.ndarray
        The interpolated scalar value as a zero-dimensional array.

    Examples
    --------
    >>> import numpy as np
    >>> from psi_io.psi_io import interpolate_point_from_1d_slice
    >>> x = np.linspace(0.0, 1.0, 11)
    >>> f = x ** 2
    >>> float(interpolate_point_from_1d_slice(0.5, x, f))
    0.25
    """
    if scalex[0] > scalex[-1]:
        scalex, values = scalex[::-1], values[::-1]
    xi_ = int(np.searchsorted(scalex, xi))
    sx_ = slice(*_check_index_ranges(len(scalex), xi_, xi_))
    return _np_linear_interpolation([xi], [scalex[sx_]], values[sx_])


def interpolate_point_from_2d_slice(xi, yi, scalex, scaley, values):
    """
    Interpolate a point from a 2D slice using bilinear interpolation.

    Parameters
    ----------
    xi : float
        The coordinate value for the first dimension.
    yi : float
        The coordinate value for the second dimension.
    scalex : np.ndarray
        The coordinate scale array for the first dimension.
    scaley : np.ndarray
        The coordinate scale array for the second dimension.
    values : np.ndarray
        The 2D data array to interpolate, with shape ``(nx, ny)``.

    Returns
    -------
    out : np.ndarray
        The interpolated scalar value as a zero-dimensional array.

    Examples
    --------
    >>> import numpy as np
    >>> from psi_io.psi_io import interpolate_point_from_2d_slice
    >>> x = np.linspace(0.0, 1.0, 11)
    >>> y = np.linspace(0.0, 1.0, 11)
    >>> f = np.outer(x, y)
    >>> float(interpolate_point_from_2d_slice(0.5, 0.5, x, y, f))
    0.25
    """
    values = np.transpose(values)
    if scalex[0] > scalex[-1]:
        scalex, values = scalex[::-1], values[::-1, :]
    if scaley[0] > scaley[-1]:
        scaley, values = scaley[::-1], values[:, ::-1]
    xi_, yi_ = int(np.searchsorted(scalex, xi)), int(np.searchsorted(scaley, yi))
    sx_, sy_ = slice(*_check_index_ranges(len(scalex), xi_, xi_)), slice(*_check_index_ranges(len(scaley), yi_, yi_))
    return _np_bilinear_interpolation([xi, yi], [scalex[sx_], scaley[sy_]], values[(sx_, sy_)])


def _rdhdf_nd(hdf_filename: str,
              dimensionality: int
              ) -> Tuple[np.ndarray, ...]:
    """Read an n-dimensional PSI-style HDF file; shared implementation for `rdhdf_1d/2d/3d`.

    Parameters
    ----------
    hdf_filename : PathLike
        The path to the HDF5 (.h5) or HDF4 (.hdf) file to read.
    dimensionality : int
        The expected number of dimensions of the dataset.

    Returns
    -------
    out : tuple[np.ndarray, ...]
        The coordinate scales for each dimension followed by the data array.

    Raises
    ------
    ValueError
        If the dataset dimensionality does not match ``dimensionality``.

    Examples
    --------
    >>> from psi_io.psi_io import _rdhdf_nd
    >>> from psi_data import fetch_mas_data
    >>> r, t, p, f = _rdhdf_nd(fetch_mas_data().cor_br, dimensionality=3)
    >>> f.shape
    (299, 142, 255)
    """
    f, *scales = read_hdf_data(hdf_filename)
    if f.ndim != dimensionality:
        err = f'Expected {dimensionality}D data, got {f.ndim}D data instead.'
        raise ValueError(err)
    scales = scales or (np.empty(0) for _ in f.shape)
    return *scales, f


def _wrhdf_nd(hdf_filename: str,
              data: np.ndarray,
              *scales: Sequence[Union[np.ndarray, None]],
              dimensionality: int,
              sync_dtype: bool = True,
              **kwargs
              ) -> None:
    """Write an n-dimensional PSI-style HDF file; shared implementation for `wrhdf_1d/2d/3d`.

    Parameters
    ----------
    hdf_filename : PathLike
        The path to the HDF5 (.h5) or HDF4 (.hdf) file to write.
    data : np.ndarray
        The data array to write.
    *scales : Sequence[np.ndarray | None]
        The scales (coordinate arrays) for each dimension.
    dimensionality : int
        The expected number of dimensions of the data array.
    sync_dtype : bool
        If True, scale dtypes are cast to match the data array dtype.
    **kwargs
        Additional keyword arguments forwarded to :func:`write_hdf_data`.

    Raises
    ------
    ValueError
        If the data dimensionality does not match ``dimensionality``.

    Examples
    --------
    >>> import tempfile, numpy as np
    >>> from pathlib import Path
    >>> from psi_io.psi_io import _wrhdf_nd
    >>> from psi_io import read_hdf_data
    >>> f = np.ones((10, 20, 30), dtype=np.float32)
    >>> r = np.linspace(1.0, 5.0, 10, dtype=np.float32)
    >>> t = np.linspace(0.0, 3.14, 20, dtype=np.float32)
    >>> p = np.linspace(0.0, 6.28, 30, dtype=np.float32)
    >>> with tempfile.TemporaryDirectory() as d:
    ...     _wrhdf_nd(Path(d) / "out.h5", f, r, t, p, dimensionality=3)
    ...     data, *_ = read_hdf_data(Path(d) / "out.h5")
    ...     data.shape
    (10, 20, 30)
    """
    if data.ndim != dimensionality:
        err = f'Expected {dimensionality}D data, got {data.ndim}D data instead.'
        raise ValueError(err)
    write_hdf_data(hdf_filename, data, *scales, sync_dtype=sync_dtype, **kwargs)


def _get_scales_nd_h5(ifile: Union[ Path, str], /,
                      dimensionality: int,
                      dataset_id: Optional[str] = None,
                      ):
    """HDF5 (.h5) version of :func:`get_scales_1d` / :func:`get_scales_2d` / :func:`get_scales_3d`.

    Parameters
    ----------
    ifile : PathLike
        The path to the HDF5 (.h5) file.
    dimensionality : int
        The expected number of dimensions of the dataset.
    dataset_id : str, optional
        The identifier of the dataset whose scales are read.
        If ``None``, the default PSI dataset (``'Data'``) is used.

    Returns
    -------
    out : tuple[np.ndarray, ...]
        The coordinate scale arrays for each dimension.

    Raises
    ------
    ValueError
        If the dataset dimensionality does not match ``dimensionality``.
    ValueError
        If any dimension has no associated scale.

    Examples
    --------
    >>> from psi_io.psi_io import _get_scales_nd_h5
    >>> from psi_data import fetch_mas_data
    >>> r, t, p = _get_scales_nd_h5(fetch_mas_data().cor_br, dimensionality=3)
    >>> r.shape, t.shape, p.shape
    ((255,), (142,), (299,))
    """
    with h5.File(ifile, 'r') as hdf:
        data = hdf[dataset_id or PSI_DATA_ID['h5']]
        ndim = data.ndim
        if ndim != dimensionality:
            err = f'Expected {dimensionality}D data, got {ndim}D data instead.'
            raise ValueError(err)
        scales = []
        for dim in data.dims:
            if dim:
                scales.append(dim[0][:])
            else:
                raise ValueError(f'Dimension has no scale associated with it.')
    return tuple(scales)


def _get_scales_nd_h4(ifile: Union[ Path, str], /,
                      dimensionality: int,
                      dataset_id: Optional[str] = None,
                      ):
    """HDF4 (.hdf) version of :func:`get_scales_1d` / :func:`get_scales_2d` / :func:`get_scales_3d`.

    Parameters
    ----------
    ifile : PathLike
        The path to the HDF4 (.hdf) file.
    dimensionality : int
        The expected number of dimensions of the dataset.
    dataset_id : str, optional
        The identifier of the dataset whose scales are read.
        If ``None``, the default PSI dataset (``'Data-Set-2'``) is used.

    Returns
    -------
    out : tuple[np.ndarray, ...]
        The coordinate scale arrays for each dimension.

    Raises
    ------
    ValueError
        If the dataset dimensionality does not match ``dimensionality``.
    ValueError
        If any dimension has no associated scale.

    Examples
    --------
    >>> from psi_io.psi_io import _get_scales_nd_h4
    >>> from psi_data import fetch_mas_data
    >>> path = fetch_mas_data(hdf=4).cor_br  # doctest: +SKIP
    >>> r, t, p = _get_scales_nd_h4(path, dimensionality=3)  # doctest: +SKIP
    >>> r.shape  # doctest: +SKIP
    (151,)
    """
    hdf = h4.SD(str(ifile))
    data = hdf.select(dataset_id or PSI_DATA_ID['h4'])
    ndim = data.info()[1]
    if ndim != dimensionality:
        err = f'Expected {dimensionality}D data, got {ndim}D data instead.'
        raise ValueError(err)
    scales = []
    for k_, v_ in reversed(data.dimensions(full=1).items()):
        if v_[3]:
            scales.append(hdf.select(k_)[:])
        else:
            raise ValueError('Dimension has no scale associated with it.')
    return tuple(scales)


def _read_h5_meta(ifile: PathLike, /,
                  dataset_id: Optional[str] = None
                  ):
    """HDF5 (.h5) version of :func:`read_hdf_meta`.

    Examples
    --------
    >>> from psi_io.psi_io import _read_h5_meta
    >>> from psi_data import fetch_mas_data
    >>> meta = _read_h5_meta(fetch_mas_data().cor_br)
    >>> meta[0].name
    'Data'
    """
    with h5.File(ifile, 'r') as hdf:
        # Raises KeyError if ``dataset_id`` not found
        # If ``dataset_id`` is None, get all non-scale :class:`h5.Dataset`s
        if dataset_id:
            datasets = (dataset_id, hdf[dataset_id]),
        else:
            datasets = ((k, v) for k, v in hdf.items() if not v.is_scale)

        # One should avoid multiple calls to ``dimproxy[0]`` – *e.g.* ``dimproxy[0].dtype`` and
        # ``dimproxy[0].shape`` – because the __getitem__ method creates and returns a new
        # :class:`~h5.DimensionProxy` object each time it is called. [Does this matter? Probably not.]
        return [HdfDataMeta(name=k,
                            type=v.dtype,
                            shape=v.shape,
                            # Exclude HDF5 dimension-scale bookkeeping attributes
                            # (``DIMENSION_LIST`` / ``DIMENSION_LABELS``).  These are
                            # managed automatically by ``attach_scale`` and hold object
                            # references that are not real user attributes; surfacing
                            # them breaks round-trips such as :func:`convert` to HDF4.
                            attr={ak: av for ak, av in v.attrs.items()
                                  if not ak.startswith("DIMENSION")},
                            scales=[HdfScaleMeta(name=dimproxy.label,
                                                 type=dim.dtype,
                                                 shape=dim.shape,
                                                 attr=dict(dim.attrs),
                                                 imin=dim[0],
                                                 imax=dim[-1])
                                    for dimproxy in v.dims if dimproxy and (dim := dimproxy[0])])
                for k, v in datasets]


def _read_h4_meta(ifile: PathLike, /,
                  dataset_id: Optional[str] = None
                  ):
    """HDF4 (.hdf) version of :func:`read_hdf_meta`.

    Examples
    --------
    >>> from psi_io.psi_io import _read_h4_meta
    >>> from psi_data import fetch_mas_data
    >>> meta = _read_h4_meta(fetch_mas_data(hdf=4).cor_br)  # doctest: +SKIP
    >>> meta[0].name  # doctest: +SKIP
    'Data-Set-2'
    """
    hdf = h4.SD(str(ifile))
    # Raises HDF4Error if ``dataset_id`` not found
    # If ``dataset_id`` is None, get all non-scale :class:`pyhdf.SD.SDS`s
    if dataset_id:
        datasets = (dataset_id, hdf.select(dataset_id)),
    else:
        datasets = ((k, hdf.select(k)) for k in hdf.datasets().keys() if not hdf.select(k).iscoordvar())

    # The inner list comprehension differs in approach from the HDF5 version because calling
    # ``dimensions(full=1)`` on an :class:`~pyhdf.SD.SDS` returns a dictionary of dimension
    # dataset identifiers (keys) and tuples containing dimension metadata (values). Even if no
    # coordinate-variable datasets are defined, this dictionary is still returned; the only
    # indication that the datasets returned do not exist is that the "type" field (within the
    # tuple of dimension metadata) is set to 0.

    # Also, one cannot avoid multiple calls to ``hdf.select(k_)`` within the inner list comprehension
    # because :class:`~pyhdf.SD.SDS` objects do not define a ``__bool__`` method, and the fallback
    # behavior of Python is to assess if the __len__ method returns a non-zero value (which, in
    # this case, always returns 0).
    return [HdfDataMeta(name=k,
                        type=SDC_TYPE_CONVERSIONS[v.info()[3]],
                        shape=_cast_shape_tuple(v.info()[2]),
                        attr=v.attributes(),
                        scales=[HdfScaleMeta(name=k_,
                                             type=SDC_TYPE_CONVERSIONS[v_[3]],
                                             shape=_cast_shape_tuple(v_[0]),
                                             attr=hdf.select(k_).attributes(),
                                             imin=hdf.select(k_)[0],
                                             imax=hdf.select(k_)[-1])
                                for k_, v_ in v.dimensions(full=1).items() if v_[3]])
            for k, v in datasets]


def _read_h5_rtp(ifile: Union[ Path, str], /):
    """HDF5 (.h5) version of :func:`read_rtp_meta`.

    Examples
    --------
    >>> from psi_io.psi_io import _read_h5_rtp
    >>> from psi_data import fetch_mas_data
    >>> meta = _read_h5_rtp(fetch_mas_data().cor_br)
    >>> sorted(meta.keys())
    ['p', 'r', 't']
    """
    with h5.File(ifile, 'r') as hdf:
        return {k: (hdf[v].size, hdf[v][0], hdf[v][-1])
                for k, v in zip('rtp', PSI_SCALE_ID['h5'])}


def _read_h4_rtp(ifile: Union[ Path, str], /):
    """HDF4 (.hdf) version of :func:`read_rtp_meta`.

    Examples
    --------
    >>> from psi_io.psi_io import _read_h4_rtp
    >>> from psi_data import fetch_mas_data
    >>> meta = _read_h4_rtp(fetch_mas_data(hdf=4).cor_br)  # doctest: +SKIP
    >>> sorted(meta.keys())  # doctest: +SKIP
    ['p', 'r', 't']
    """
    hdf = h4.SD(str(ifile))
    return {k: (hdf.select(v).info()[2], hdf.select(v)[0], hdf.select(v)[-1])
            for k, v in zip('ptr', PSI_SCALE_ID['h4'])}


def _read_h5_data(ifile: PathLike, /,
                  dataset_id: Optional[str] = None,
                  return_scales: bool = True,
                  ) -> Union[np.ndarray, Tuple[np.ndarray]]:
    """HDF5 (.h5) version of :func:`read_hdf_data`.

    Examples
    --------
    >>> from psi_io.psi_io import _read_h5_data
    >>> from psi_data import fetch_mas_data
    >>> data, r, t, p = _read_h5_data(fetch_mas_data().cor_br)
    >>> data.shape
    (299, 142, 255)
    """
    with h5.File(ifile, 'r') as hdf:
        data = hdf[dataset_id or PSI_DATA_ID['h5']]
        dataset = data[:]
        if return_scales:
            scales = [dim[0][:] for dim in data.dims if dim]
            return dataset, *scales
        return dataset


def _read_h4_data(ifile: PathLike, /,
                  dataset_id: Optional[str] = None,
                  return_scales: bool = True,
                  ) -> Union[np.ndarray, Tuple[np.ndarray]]:
    """HDF4 (.hdf) version of :func:`read_hdf_data`.

    Examples
    --------
    >>> from psi_io.psi_io import _read_h4_data
    >>> from psi_data import fetch_mas_data
    >>> data, *_ = _read_h4_data(fetch_mas_data(hdf=4).cor_br)  # doctest: +SKIP
    >>> data.shape  # doctest: +SKIP
    (299, 142, 255)
    """
    hdf = h4.SD(str(ifile))
    data = hdf.select(dataset_id or PSI_DATA_ID['h4'])
    if return_scales:
        out = (data[:],
               *[hdf.select(k_)[:] for k_, v_ in reversed(data.dimensions(full=1).items()) if v_[3]])
    else:
        out = data[:]
    return out


def _read_h5_by_index(ifile: PathLike, /,
                      *xi: Union[int, Tuple[Union[int, None], Union[int, None]], None],
                      dataset_id: Optional[str] = None,
                      return_scales: bool = True,
                      ) -> Union[np.ndarray, Tuple[np.ndarray]]:
    """HDF5 (.h5) version of :func:`read_hdf_by_index`.

    Examples
    --------
    >>> from psi_io.psi_io import _read_h5_by_index
    >>> from psi_data import fetch_mas_data
    >>> f, r, t, p = _read_h5_by_index(fetch_mas_data().cor_br, 0, None, None)
    >>> f.shape
    (299, 142, 1)
    """
    with h5.File(ifile, 'r') as hdf:
        data = hdf[dataset_id or PSI_DATA_ID['h5']]
        if len(xi) != data.ndim:
            raise ValueError(f"len(xi) must equal the number of scales for {dataset_id}")
        slices = [_parse_index_inputs(slice_input) for slice_input in xi]
        dataset = data[tuple(reversed(slices))]
        if return_scales:
            scales = [dim[0][si] for si, dim in zip(slices, data.dims) if dim]
            return dataset, *scales
        return dataset

def _read_h4_by_index(ifile: PathLike, /,
                      *xi: Union[int, Tuple[Union[int, None], Union[int, None]], None],
                      dataset_id: Optional[str] = None,
                      return_scales: bool = True,
                      ) -> Union[np.ndarray, Tuple[np.ndarray]]:
    """HDF4 (.hdf) version of :func:`read_hdf_by_index`.

    Examples
    --------
    >>> from psi_io.psi_io import _read_h4_by_index
    >>> from psi_data import fetch_mas_data
    >>> path = fetch_mas_data(hdf=4).cor_br  # doctest: +SKIP
    >>> f, r, t, p = _read_h4_by_index(path, 0, None, None)  # doctest: +SKIP
    >>> f.shape  # doctest: +SKIP
    (299, 142, 1)
    """
    hdf = h4.SD(str(ifile))
    data = hdf.select(dataset_id or PSI_DATA_ID['h4'])
    ndim = data.info()[1]
    if len(xi) != ndim:
        raise ValueError(f"len(xi) must equal the number of scales for {dataset_id}")
    slices = [_parse_index_inputs(slice_input) for slice_input in xi]
    dataset = data[tuple(reversed(slices))]
    if return_scales:
        scales = [hdf.select(k_)[si] for si, (k_, v_) in zip(slices, reversed(data.dimensions(full=1).items())) if v_[3]]
        return dataset, *scales
    return dataset


def _read_h5_by_value(ifile: PathLike, /,
                      *xi: Union[float, Tuple[float, float], None],
                      dataset_id: Optional[str] = None,
                      return_scales: bool = True,
                      ) -> Union[np.ndarray, Tuple[np.ndarray]]:
    """HDF5 (.h5) version of :func:`read_hdf_by_value`.

    Examples
    --------
    >>> from psi_io.psi_io import _read_h5_by_value
    >>> from psi_data import fetch_mas_data
    >>> f, r, t, p = _read_h5_by_value(fetch_mas_data().cor_br, 15, None, None)
    >>> f.shape
    (299, 142, 2)
    """
    with h5.File(ifile, 'r') as hdf:
        data = hdf[dataset_id or PSI_DATA_ID['h5']]
        if len(xi) != data.ndim:
            raise ValueError(f"len(xi) must equal the number of scales for {dataset_id}")
        slices = []
        for dimproxy, value in zip(data.dims, xi):
            if dimproxy:
                slices.append(_parse_value_inputs(dimproxy[0], value))
            elif value is None:
                slices.append(slice(None))
            else:
                raise ValueError("Cannot slice by value on dimension without scales")
        dataset = data[tuple(reversed(slices))]
        if return_scales:
            scales = [dim[0][si] for si, dim in zip(slices, data.dims) if dim]
            return dataset, *scales
        return dataset


def _read_h4_by_value(ifile: PathLike, /,
                      *xi: Union[float, Tuple[float, float], None],
                      dataset_id: Optional[str] = None,
                      return_scales: bool = True,
                      ) -> Union[np.ndarray, Tuple[np.ndarray]]:
    """HDF4 (.hdf) version of :func:`read_hdf_by_value`.

    Examples
    --------
    >>> from psi_io.psi_io import _read_h4_by_value
    >>> from psi_data import fetch_mas_data
    >>> path = fetch_mas_data(hdf=4).cor_br  # doctest: +SKIP
    >>> f, r, t, p = _read_h4_by_value(path, 15, None, None)  # doctest: +SKIP
    >>> f.shape  # doctest: +SKIP
    (299, 142, 2)
    """
    hdf = h4.SD(str(ifile))
    data = hdf.select(dataset_id or PSI_DATA_ID['h4'])
    ndim = data.info()[1]
    if len(xi) != ndim:
        raise ValueError(f"len(xi) must equal the number of scales for {dataset_id}")
    slices = []
    for (k_, v_), value in zip(reversed(data.dimensions(full=1).items()), xi):
        if v_[3] != 0:
            slices.append(_parse_value_inputs(hdf.select(k_), value))
        elif value is None:
            slices.append(slice(None))
        else:
            raise ValueError("Cannot slice by value on dimension without scales")
    dataset = data[tuple(reversed(slices))]
    if return_scales:
        scales = [hdf.select(k_)[si] for si, (k_, v_) in zip(slices, reversed(data.dimensions(full=1).items())) if v_[3]]
        return dataset, *scales
    return dataset


def _read_h5_by_ivalue(ifile: PathLike, /,
                       *xi: Union[float, Tuple[float, float], None],
                       dataset_id: Optional[str] = None,
                       return_scales: bool = True,
                       ) -> Union[np.ndarray, Tuple[np.ndarray]]:
    """HDF5 (.h5) version of :func:`read_hdf_by_ivalue`.

    Examples
    --------
    >>> from psi_io.psi_io import _read_h5_by_ivalue
    >>> from psi_data import fetch_mas_data
    >>> f, r_idx, t_idx, p_idx = _read_h5_by_ivalue(fetch_mas_data().cor_br, 2.7, None, None)
    >>> r_idx.shape
    (2,)
    """
    with h5.File(ifile, 'r') as hdf:
        data = hdf[dataset_id or PSI_DATA_ID['h5']]
        if len(xi) != data.ndim:
            raise ValueError(f"len(xi) must equal the number of scales for {dataset_id}")
        slices = [_parse_ivalue_inputs(*args) for args in zip(reversed(data.shape), xi)]
        dataset = data[tuple(reversed(slices))]
        if return_scales:
            scales = [np.arange(si.start or 0, si.stop or size) for si, size in zip(slices, reversed(data.shape))]
            return dataset, *scales
        return dataset


def _read_h4_by_ivalue(ifile: PathLike, /,
                       *xi: Union[float, Tuple[float, float], None],
                       dataset_id: Optional[str] = None,
                       return_scales: bool = True,
                       ) -> Union[np.ndarray, Tuple[np.ndarray]]:
    """HDF4 (.hdf) version of :func:`read_hdf_by_ivalue`.

    Examples
    --------
    >>> from psi_io.psi_io import _read_h4_by_ivalue
    >>> from psi_data import fetch_mas_data
    >>> path = fetch_mas_data(hdf=4).cor_br  # doctest: +SKIP
    >>> f, r_idx, t_idx, p_idx = _read_h4_by_ivalue(path, 2.7, None, None)  # doctest: +SKIP
    >>> r_idx.shape  # doctest: +SKIP
    (2,)
    """
    hdf = h4.SD(str(ifile))
    data = hdf.select(dataset_id or PSI_DATA_ID['h4'])
    ndim, shape = data.info()[1], _cast_shape_tuple(data.info()[2])
    if len(xi) != ndim:
        raise ValueError(f"len(xi) must equal the number of scales for {dataset_id}")
    slices = [_parse_ivalue_inputs(*args) for args in zip(reversed(shape), xi)]
    dataset = data[tuple(reversed(slices))]
    if return_scales:
        scales = [np.arange(si.start or 0, si.stop or size) for si, size in zip(slices, reversed(shape))]
        return dataset, *scales
    return dataset


def _write_h4_data(ifile: PathLike, /,
                   data: np.ndarray,
                   *scales: Sequence[np.ndarray],
                   dataset_id: Optional[str] = None,
                   sync_dtype: bool = False,
                   strict: bool = True,
                   **kwargs) -> Path:
    """HDF4 (.hdf) version of :func:`write_hdf_data`.

    Examples
    --------
    >>> import tempfile, numpy as np
    >>> from pathlib import Path
    >>> from psi_io.psi_io import _write_h4_data, _read_h4_data
    >>> f = np.ones((10,), dtype=np.float32)
    >>> x = np.linspace(0.0, 1.0, 10, dtype=np.float32)
    >>> with tempfile.TemporaryDirectory() as d:  # doctest: +SKIP
    ...     _write_h4_data(Path(d) / "out.hdf", f, x)  # doctest: +SKIP
    ...     data, *_ = _read_h4_data(Path(d) / "out.hdf")  # doctest: +SKIP
    ...     data.shape  # doctest: +SKIP
    (10,)
    """
    dataid = dataset_id or PSI_DATA_ID['h4']
    h4file = h4.SD(str(ifile), h4.SDC.WRITE | h4.SDC.CREATE | h4.SDC.TRUNC)
    sds_id = h4file.create(dataid, _dtype_to_sdc(data.dtype), data.shape)

    if scales:
        for i, scale in enumerate(reversed(scales)):
            if scale is not None:
                if sync_dtype:
                    scale = scale.astype(data.dtype)
                sds_id.dim(i).setscale(_dtype_to_sdc(scale.dtype), scale.tolist())

    if kwargs:
        for k, v in kwargs.items():
            npv = np.asarray(v)
            attr_ = sds_id.attr(k)
            try:
                val = npv.tolist()
                if isinstance(val, bytes):
                    val = val.decode('latin-1')
                attr_.set(_dtype_to_sdc(npv.dtype), val)
            except KeyError as e:
                if strict:
                    raise KeyError(f"Failed to set attribute '{k}' on dataset '{dataid}'") from e
                else:
                    print(f"Warning: Failed to set attribute '{k}' on dataset '{dataid}'; skipping.")

    sds_id.set(data)
    sds_id.endaccess()
    h4file.end()

    return ifile


def _write_h5_data(ifile: PathLike, /,
                   data: np.ndarray,
                   *scales: Sequence[np.ndarray],
                   dataset_id: Optional[str] = None,
                   sync_dtype: bool = False,
                   strict: bool = True,
                   **kwargs) -> Path:
    """HDF5 (.h5) version of :func:`write_hdf_data`.

    Examples
    --------
    >>> import tempfile, numpy as np
    >>> from pathlib import Path
    >>> from psi_io.psi_io import _write_h5_data, _read_h5_data
    >>> f = np.ones((10,), dtype=np.float32)
    >>> x = np.linspace(0.0, 1.0, 10, dtype=np.float32)
    >>> with tempfile.TemporaryDirectory() as d:
    ...     _ = _write_h5_data(Path(d) / "out.h5", f, x)
    ...     data, *_ = _read_h5_data(Path(d) / "out.h5")
    ...     data.shape
    (10,)
    """
    dataid = dataset_id or PSI_DATA_ID['h5']
    with h5.File(ifile, "w") as h5file:
        dataset = h5file.create_dataset(dataid, data=data, dtype=data.dtype, shape=data.shape)

        if scales:
            for i, scale in enumerate(scales):
                if scale is not None:
                    if sync_dtype:
                        scale = scale.astype(data.dtype)
                    h5file.create_dataset(f"dim{i+1}", data=scale, dtype=scale.dtype, shape=scale.shape)
                    h5file[dataid].dims[i].attach_scale(h5file[f"dim{i+1}"])
                    h5file[dataid].dims[i].label = f"dim{i+1}"

        if kwargs:
            for key, value in kwargs.items():
                if key.startswith('DIMENSION'):
                    # Skip HDF5 dimension-scale bookkeeping attributes —
                    # these are managed by attach_scale above and must not
                    # be overwritten with stale object references.
                    continue
                try:
                    dataset.attrs[key] = value
                except TypeError as e:
                    if strict:
                        raise TypeError(f"Failed to set attribute '{key}' on dataset '{dataid}'") from e
                    else:
                        print(f"Warning: Failed to set attribute '{key}' on dataset '{dataid}'; skipping.")

    return ifile


def _write_h4_meta(ifile: PathLike, /,
                   meta: Optional[Mapping[str, Mapping[str, Any]]] = None,
                   **kwargs) -> Path:
    """HDF4 (.hdf) version of :func:`write_hdf_meta`."""
    metadata = dict(meta or {})
    h4file = h4.SD(str(ifile), h4.SDC.READ | h4.SDC.WRITE)

    for k, v in kwargs.items():
        npv = np.asarray(v)
        attr_ = h4file.attr(k)
        val = npv.tolist()
        if isinstance(val, bytes):
            val = val.decode('latin-1')
        attr_.set(_dtype_to_sdc(npv.dtype), val)
    for key, value in metadata.items():
        sds_id = h4file.select(key)
        for k, v in dict(value).items():
            npv = np.asarray(v)
            attr_ = sds_id.attr(k)
            val = npv.tolist()
            if isinstance(val, bytes):
                val = val.decode('latin-1')
            attr_.set(_dtype_to_sdc(npv.dtype), val)
    h4file.end()

    return ifile


def _write_h5_meta(ifile: PathLike, /,
                   meta: Optional[Mapping[str, Mapping[str, Any]]] = None,
                   **kwargs) -> Path:
    """HDF5 (.h5) version of :func:`write_hdf_meta`."""
    metadata = dict(meta or {})
    with h5.File(ifile, "r+") as h5file:
        if kwargs:
            h5file.attrs.update(**kwargs)
        for key, value in metadata.items():
            h5file[key].attrs.update(**dict(value))
    return ifile


def _np_linear_interpolation(xi: Sequence, scales: Sequence, values: np.ndarray):
    """
    Perform linear interpolation over one dimension.

    Parameters
    ----------
    xi : list[float | None]
        Target coordinate values for each dimension; ``None`` marks a
        "free" (not interpolated) dimension.
    scales : list[np.ndarray]
        Coordinate scale arrays, one per dimension.
    values : np.ndarray
        The data array to interpolate.

    Returns
    -------
    out : np.ndarray
        The interpolated data array.

    Examples
    --------
    >>> import numpy as np
    >>> from psi_io.psi_io import _np_linear_interpolation
    >>> scales = [np.array([0.0, 1.0])]
    >>> values = np.array([0.0, 2.0])
    >>> float(_np_linear_interpolation([0.5], scales, values))
    1.0
    """
    index0 = next((i for i, v in enumerate(xi) if v is not None), None)
    t = (xi[index0] - scales[index0][0])/(scales[index0][1] - scales[index0][0])
    f0 = [slice(None, None)]*values.ndim
    f1 = [slice(None, None)]*values.ndim
    f0[index0] = 0
    f1[index0] = 1

    return (1 - t)*values[tuple(f0)] + t*values[tuple(f1)]


def _np_bilinear_interpolation(xi, scales, values):
    """
    Perform bilinear interpolation over two dimensions.

    Parameters
    ----------
    xi : Sequence[float | None]
        Target coordinate values for each dimension; ``None`` marks a
        "free" (not interpolated) dimension.
    scales : Sequence[np.ndarray]
        Coordinate scale arrays, one per dimension.
    values : np.ndarray
        The data array to interpolate.

    Returns
    -------
    out : np.ndarray
        The interpolated data array.

    Examples
    --------
    >>> import numpy as np
    >>> from psi_io.psi_io import _np_bilinear_interpolation
    >>> scales = [np.array([0.0, 1.0]), np.array([0.0, 1.0])]
    >>> values = np.array([[0.0, 0.0], [1.0, 2.0]])
    >>> float(_np_bilinear_interpolation([0.5, 0.5], scales, values))
    0.75
    """
    index0, index1 = [i for i, v in enumerate(xi) if v is not None]
    t, u = [(xi[i] - scales[i][0])/(scales[i][1] - scales[i][0]) for i in (index0, index1)]

    f00 = [slice(None, None)]*values.ndim
    f10 = [slice(None, None)]*values.ndim
    f01 = [slice(None, None)]*values.ndim
    f11 = [slice(None, None)]*values.ndim
    f00[index0], f00[index1] = 0, 0
    f10[index0], f10[index1] = 1, 0
    f01[index0], f01[index1] = 0, 1
    f11[index0], f11[index1] = 1, 1

    return (
          (1 - t)*(1 - u)*values[tuple(f00)] +
          t*(1 - u)*values[tuple(f10)] +
          (1 - t)*u*values[tuple(f01)] +
          t*u*values[tuple(f11)]
    )


def _np_trilinear_interpolation(xi, scales, values):
    """
    Perform trilinear interpolation over three dimensions.

    Parameters
    ----------
    xi : list[float | None]
        Target coordinate values for each dimension; ``None`` marks a
        "free" (not interpolated) dimension.
    scales : list[np.ndarray]
        Coordinate scale arrays, one per dimension.
    values : np.ndarray
        The data array to interpolate.

    Returns
    -------
    out : np.ndarray
        The interpolated data array.

    Examples
    --------
    >>> import numpy as np
    >>> from psi_io.psi_io import _np_trilinear_interpolation
    >>> scales = [np.array([0.0, 1.0]), np.array([0.0, 1.0]), np.array([0.0, 1.0])]
    >>> values = np.zeros((2, 2, 2))
    >>> values[1, 1, 1] = 8.0
    >>> float(_np_trilinear_interpolation([0.5, 0.5, 0.5], scales, values))
    1.0
    """
    index0, index1, index2 = [i for i, v in enumerate(xi) if v is not None]
    t, u, v = [(xi[i] - scales[i][0])/(scales[i][1] - scales[i][0]) for i in (index0, index1, index2)]

    f000 = [slice(None, None)]*values.ndim
    f100 = [slice(None, None)]*values.ndim
    f010 = [slice(None, None)]*values.ndim
    f110 = [slice(None, None)]*values.ndim
    f001 = [slice(None, None)]*values.ndim
    f101 = [slice(None, None)]*values.ndim
    f011 = [slice(None, None)]*values.ndim
    f111 = [slice(None, None)]*values.ndim

    f000[index0], f000[index1], f000[index2] = 0, 0, 0
    f100[index0], f100[index1], f100[index2] = 1, 0, 0
    f010[index0], f010[index1], f010[index2] = 0, 1, 0
    f110[index0], f110[index1], f110[index2] = 1, 1, 0
    f001[index0], f001[index1], f001[index2] = 0, 0, 1
    f101[index0], f101[index1], f101[index2] = 1, 0, 1
    f011[index0], f011[index1], f011[index2] = 0, 1, 1
    f111[index0], f111[index1], f111[index2] = 1, 1, 1

    c00 = values[tuple(f000)]*(1 - t) + values[tuple(f100)]*t
    c10 = values[tuple(f010)]*(1 - t) + values[tuple(f110)]*t
    c01 = values[tuple(f001)]*(1 - t) + values[tuple(f101)]*t
    c11 = values[tuple(f011)]*(1 - t) + values[tuple(f111)]*t

    c0 = c00*(1 - u) + c10*u
    c1 = c01*(1 - u) + c11*u

    return c0*(1 - v) + c1*v


def _check_index_ranges(arr_size: int,
                        i0: Union[int, np.integer],
                        i1: Union[int, np.integer]
                        ) -> Tuple[int, int]:
    """
    Adjust index ranges to ensure they cover at least two indices.

    Parameters
    ----------
    arr_size : int
        The size of the array along the dimension.
    i0 : int
        The starting index.
    i1 : int
        The ending index.

    Returns
    -------
    out : tuple[int, int]
        Adjusted starting and ending indices.

    Notes
    -----
    This function ensures that the range between `i0` and `i1` includes at least
    two indices for interpolation purposes.

    Examples
    --------
    >>> from psi_io.psi_io import _check_index_ranges
    >>> _check_index_ranges(10, 0, 0)
    (0, 2)
    >>> _check_index_ranges(10, 5, 5)
    (4, 6)
    >>> _check_index_ranges(10, 10, 10)
    (8, 10)
    """
    i0, i1 = int(i0), int(i1)
    if i0 == 0:
        return (i0, i1 + 2) if i1 == 0 else (i0, i1 + 1)
    elif i0 == arr_size:
        return i0 - 2, i1
    else:
        return i0 - 1, i1 + 1


def _cast_shape_tuple(input: Union[int, Sequence[int]]
                      ) -> tuple[int, ...]:
    """
    Cast an input to a tuple of integers.

    Parameters
    ----------
    input : int | Sequence[int]
        The input to cast.

    Returns
    -------
    out : tuple[int, ...]
        The input cast as a tuple of integers.

    Raises
    ------
    TypeError
        If the input is neither an integer nor an iterable of integers.

    Examples
    --------
    >>> from psi_io.psi_io import _cast_shape_tuple
    >>> _cast_shape_tuple(5)
    (5,)
    >>> _cast_shape_tuple([3, 4, 5])
    (3, 4, 5)
    """
    if isinstance(input, int):
        return (input,)
    elif isinstance(input, Sequence):
        return tuple(int(i) for i in input)
    else:
        raise TypeError("Input must be an integer or an iterable of integers.")


def _parse_index_inputs(input: Union[int, slice, Sequence[Union[int, None]], None]
                        ) -> slice:
    """
    Parse various slice input formats into a standard slice object.

    Parameters
    ----------
    input : int | slice | Sequence[int | None] | None
        The index specification to parse:

        - *int* — selects a single-element slice ``[n, n+1)``.
        - *Sequence[int | None]* — unpacked directly into :class:`slice`.
        - ``None`` — selects the entire dimension.

    Returns
    -------
    slice
        The parsed slice object.

    Raises
    ------
    TypeError
        If the input type is unsupported.

    Examples
    --------
    >>> from psi_io.psi_io import _parse_index_inputs
    >>> _parse_index_inputs(3)
    slice(3, 4, None)
    >>> _parse_index_inputs((2, 7))
    slice(2, 7, None)
    >>> _parse_index_inputs(None)
    slice(None, None, None)
    """
    if isinstance(input, int):
        return slice(input, input + 1)
    elif isinstance(input, Sequence):
        return slice(*input)
    elif input is None:
        return slice(None)
    else:
        raise TypeError("Unsupported input type for slicing.")


def _parse_value_inputs(dimproxy,
                        value,
                        scale_exists: bool = True
                        ) -> slice:
    """Parse a scale value or value range into a slice over a dimension's coordinate array.

    Parameters
    ----------
    dimproxy : array-like
        The coordinate array (scale) for the dimension, supporting ``[:]`` indexing.
    value : float | tuple[float | None, float | None] | None
        The target value or range:

        - ``None`` — select the entire dimension.
        - *float* — find the 2-element bracket ``[a, b]`` such that ``a <= value < b``.
        - *(float, float)* — find the bracket spanning the given range.
    scale_exists : bool
        Guard flag; raises :exc:`ValueError` if ``False`` and ``value`` is not ``None``.

    Returns
    -------
    slice
        A slice object that selects the appropriate indices from the coordinate array.

    Raises
    ------
    ValueError
        If ``value`` is not ``None`` and ``scale_exists`` is ``False``.

    Examples
    --------
    >>> import numpy as np
    >>> from psi_io.psi_io import _parse_value_inputs
    >>> scale = np.linspace(0.0, 10.0, 11)
    >>> _parse_value_inputs(scale, None)
    slice(None, None, None)
    >>> _parse_value_inputs(scale, 4.5)
    slice(4, 6, None)
    >>> _parse_value_inputs(scale, (2.0, 7.0))
    slice(1, 8, None)
    """
    if value is None:
        return slice(None)
    if not scale_exists:
        raise ValueError("Cannot parse value inputs when scale does not exist.")
    dim = dimproxy[:]
    if not isinstance(value, Sequence):
        insert_index = np.searchsorted(dim, value)
        return slice(*_check_index_ranges(dim.size, insert_index, insert_index))
    else:
        temp_range = list(value)
        if temp_range[0] is None:
            temp_range[0] = -np.inf
        if temp_range[-1] is None:
            temp_range[-1] = np.inf
        insert_indices = np.searchsorted(dim, temp_range)
        return slice(*_check_index_ranges(dim.size, *insert_indices))


def _parse_ivalue_inputs(dimsize,
                         input: Union[Union[int, float], slice, Sequence[Union[Union[int, float], None]], None]
                         ) -> slice:
    """
    Parse a sub-index value or range into a slice over a dimension of a given size.

    Parameters
    ----------
    dimsize : int
        The size of the array along the dimension; used to clamp the resulting slice.
    input : int | float | Sequence[int | float | None] | None
        The sub-index specification to parse:

        - ``None`` — select the entire dimension.
        - *int* or *float* (:math:`a`) — returns ``slice(floor(a), ceil(a))``
          (clamped to ``[0, dimsize - 1]``).
        - *Sequence* ``(a, b)`` — returns ``slice(floor(a), ceil(b))``
          (clamped to ``[0, dimsize - 1]``).

    Returns
    -------
    slice
        The parsed slice object, guaranteed to span at least 2 elements.

    Raises
    ------
    TypeError
        If the input type is unsupported.

    Examples
    --------
    >>> from psi_io.psi_io import _parse_ivalue_inputs
    >>> _parse_ivalue_inputs(10, None)
    slice(None, None, None)
    >>> _parse_ivalue_inputs(10, 2.7)
    slice(2, 4, None)
    >>> _parse_ivalue_inputs(10, (1.3, 4.8))
    slice(1, 5, None)
    """
    if input is None:
        return slice(None)
    elif isinstance(input, (int, float)):
        i0, i1 = math.floor(input), math.ceil(input)
    elif isinstance(input, Sequence):
        i0, i1 = math.floor(input[0]), math.ceil(input[1])
    else:
        raise TypeError("Unsupported input type for slicing.")

    if i0 > i1:
        i0, i1 = i1, i0
    i0, i1 = max(0, i0), min(dimsize - 1, i1)
    if i0 > i1:
        i0, i1 = i1, i0
    i0, i1 = max(0, i0), min(dimsize - 1, i1)
    if (i1 - i0) < 2:
        return slice(i0, i1 + 2 - (i1-i0))
    else:
        return slice(i0, i1)
