"""
runs the simulation
"""


import ctypes
from pathlib import Path
import shutil
from subprocess import run

import numpy as np

HERE = Path(__file__).parent

lib_path = HERE / "../zig-out/lib/libtf_permute.dylib"
assert lib_path.exists()


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


def simulate(
    n_cells: int, n_genes: int, prob: float, resolution: int, seed: int
) -> np.ndarray:
    # Load the shared library
    # lib_path = os.path.abspath("libtf_permute.so")
    lib = ctypes.CDLL(str(lib_path))

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
        # back to the top
        simulate(n_cells, n_genes, prob, resolution, seed)

    output_array = (ctypes.c_size_t * ctype_ncells.value)()
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


if __name__ == "__main__":
    print(simulate(118, 100, 0.2372, 10_000, 12334))
