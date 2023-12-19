const std = @import("std");
const testing = std.testing;
const ptf = @import("permute_tf.zig");
const ec = @import("exported_constants.zig");

export fn get_ncells_ngenes(n_cells: *usize, n_genes: *usize) void {
    n_cells.* = ec.EXPORTED_N_CELLS;
    n_genes.* = ec.EXPORTED_N_GENES;
}

export fn simulate_data(rng_seed: u64, prob: f64, resolution: usize, array_ptr: [*]usize) void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const alloc = gpa.allocator();
    defer {
        const deinit_status = gpa.deinit();
        if (deinit_status == .leak) @panic("Memory leak");
    }
    const result_hist = ptf.permute(ec.EXPORTED_N_GENES, resolution, ec.EXPORTED_N_CELLS, alloc, prob, rng_seed) catch unreachable;
    ptf.writeOutData(&result_hist) catch unreachable;
    for (0.., result_hist) |i, result| {
        array_ptr[i] = result;
    }
}
