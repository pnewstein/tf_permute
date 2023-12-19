"""
runs the simulation
"""

import platform
import ctypes
from pathlib import Path
import shutil
from subprocess import run
from multiprocessing import Pool
from functools import partial
from itertools import count
from ctypes import Array

import numpy as np

HERE = Path(__file__).parent


def get_lib_path() -> Path:
    path_no_suffix = HERE / "../zig-out/lib/libtf_permute"
    ext_map = {
        "darwin": ".dylib",
        "linux": ".so",
        "windows": ".dll",
    } 
    return path_no_suffix.with_suffix(ext_map[platform.system().lower()])



def recompile_zig(n_cells: int, n_genes: int):
    # write the code
    zig_code = (
        f"pub const EXPORTED_N_CELLS = {n_cells};\n"
        f"pub const EXPORTED_N_GENES = {n_genes};"
    )
    (HERE / "exported_constants.zig").write_text(zig_code, "utf-8")
    # compile the code
    zig_executable = shutil.which("zig")
    if not zig_executable:
        raise FileNotFoundError("Could not find the zig")
    print("compiling")
    run([zig_executable, "build", "-Drelease=true"], check=True, cwd=HERE)

def get_lib_output_array(n_cells: int, n_genes: int) -> tuple[ctypes.CDLL, Array[ctypes.c_size_t]]:
    try:
        lib = ctypes.CDLL(str(get_lib_path()))
    except OSError:
        recompile_zig(n_cells, n_genes)
        lib = ctypes.CDLL(str(get_lib_path()))

    # get the compiled n_genes and n_cells
    ctype_ngenes = ctypes.c_size_t()
    ctype_ncells = ctypes.c_size_t()
    get_ncells_ngenes = lib.get_ncells_ngenes
    get_ncells_ngenes.argtypes = [
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(ctypes.c_size_t),
    ]
    get_ncells_ngenes.restype = None
    get_ncells_ngenes(ctype_ncells, ctype_ngenes)

    # recompile if neccisary
    print(n_genes, ctype_ngenes.value)
    print(n_cells, ctype_ncells.value)
    if n_genes != ctype_ngenes.value or n_cells != ctype_ncells.value:
        recompile_zig(n_cells, n_genes)
        # back to the top to get lib
        return get_lib_output_array(n_cells, n_genes)
    output_array = (ctypes.c_size_t * ctype_ncells.value)()
    return lib, output_array


def simulate(
    n_cells: int, n_genes: int, prob: float, resolution: int, seed: int
) -> np.ndarray:
    # Load the shared library
    # lib_path = os.path.abspath("libtf_permute.so")
    lib, output_array = get_lib_output_array(n_cells, n_genes)

    simulate_data = lib.simulate_data

    simulate_data.argtypes = [
        ctypes.c_uint64,
        ctypes.c_double,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    simulate_data.restype = None

    # Call the function
    print("running")
    simulate_data(seed, prob, resolution, output_array)
    numpy_array = np.ctypeslib.as_array(output_array)
    return numpy_array

def simulate_until_user_input(
    n_cells: int, n_genes: int, prob: float
) -> np.ndarray:
    """
    Does simulations of resolution 100_000 until the user tells the simulation
    to stop.
    """
    _ = get_lib_output_array(n_cells, n_genes)
    chunk_size = 100_000
    simulate_partial = partial(simulate, n_cells, n_genes, prob, chunk_size)
    out = np.zeros(n_cells).astype(np.uint64)
    with Pool() as pool:
        try:
            for sim_result in pool.imap_unordered(simulate_partial, count()):
                out += sim_result
        except KeyboardInterrupt:
            # break out of the pool
            print("Keybord interupt!")
    return out


if __name__ == "__main__":
    print(simulate_until_user_input(118, 68, 0.10269192422731803))
