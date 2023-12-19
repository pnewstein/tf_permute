const FirstByte = u8;
// The ammount of spaces in each arraylist for the combo set
const FOLD_OVERALLOCATE = 3;

const std = @import("std");

fn print(a: anytype) void {
    std.debug.print("{any}\n", .{a});
}
const size = @bitSizeOf(usize);

pub fn writeOutData(data: []const usize) !void {
    // print the result
    const stdout = std.io.getStdOut();
    var buf = std.io.bufferedWriter(stdout.writer());
    var writer = buf.writer();
    for (data) |result| {
        try writer.print("{d} ", .{result});
    }
    try writer.writeByte('\n');
    try buf.flush();
}

//pub fn main() !u8 {
//var gpa = std.heap.GeneralPurposeAllocator(.{}){};
//const alloc = gpa.allocator();
//defer {
//const deinit_status = gpa.deinit();
//if (deinit_status == .leak) @panic("Memory leak");
//}
//const result_hist = try permute(10, 1_00_000, 126, alloc, 0.1, 1203);
//writeOutData(&result_hist) catch unreachable;
//return 0;
//}

pub fn permute(comptime n_genes: u16, resolution: usize, comptime n_types: u16, alloc: std.mem.Allocator, prob: f64, rng_seed: u64) ![n_types]usize {
    var fun_out = [_]usize{0} ** n_types;
    var cs = try ComboSet(n_genes).init(alloc);
    defer cs.deinit(alloc);
    var prng = std.rand.DefaultPrng.init(rng_seed);
    for (0..resolution) |_| {
        const n_intersections = try oneTrial(n_genes, n_types, alloc, prob, &cs, &prng);
        fun_out[n_intersections] += 1;
    }
    return fun_out;
}

fn ByteArray(comptime n_genes: u16) type {
    const backing_int_bits = size * try std.math.divCeil(u16, n_genes, size);
    const LastInt = @Type(.{ .Int = .{ .signedness = .unsigned, .bits = n_genes - @bitSizeOf(FirstByte) } });
    return struct {
        pub const PackedStructType = packed struct(@Type(.{ .Int = .{ .signedness = .unsigned, .bits = backing_int_bits } })) {
            firstByte: FirstByte,
            lastInt: LastInt,
            _padding: @Type(.{ .Int = .{ .signedness = .unsigned, .bits = backing_int_bits - n_genes } }) = undefined,
        };
        data: std.bit_set.ArrayBitSet(usize, n_genes),
        fn create(bool_array: [n_genes]bool) @This() {
            var data = std.bit_set.ArrayBitSet(usize, n_genes).initEmpty();
            for (0.., bool_array) |i, boolean| {
                data.setValue(i, boolean);
            }
            return @This(){
                .data = data,
            };
        }
        /// interprets the data as a byte (taking up the first 8 bits)
        fn getByte(self: @This()) FirstByte {
            const p_struct: @This().PackedStructType = @bitCast(self.data.masks);
            return p_struct.firstByte;
        }
        /// interprets the data as an integer (taking remaining bytes)
        fn getLastInt(self: @This()) LastInt {
            const p_struct: @This().PackedStructType = @bitCast(self.data.masks);
            return p_struct.lastInt;
        }
    };
}

test "create byte array" {
    const bool_array = [_]bool{ true, false, true, true, true, true, true, true };
    const ba = ByteArray(8).create(bool_array);
    for (0.., bool_array) |i, boolean| {
        try std.testing.expectEqual(ba.data.isSet(i), boolean);
    }
}

test "get first byte" {
    const ba1 = ByteArray(10).create([_]bool{false} ** 7 ++ [_]bool{true} ** 3);
    const ba2 = ByteArray(10).create([_]bool{true} ** 1 ++ [_]bool{false} ** 9);
    try std.testing.expect(ba1.getByte() != ba2.getByte());
}

test "get last int" {
    const ba1 = ByteArray(10).create([_]bool{false} ** 7 ++ [_]bool{true} ** 3);
    const ba2 = ByteArray(10).create([_]bool{true} ** 1 ++ [_]bool{false} ** 9);
    try std.testing.expect(ba1.getLastInt() != ba2.getLastInt());
}

// much like a hasmap of arraylists, but without the hash
fn ComboSet(comptime n_genes: u16) type {
    const LastInt = @Type(.{ .Int = .{ .signedness = .unsigned, .bits = n_genes - @bitSizeOf(FirstByte) } });
    const n_arrays = (std.math.maxInt(FirstByte) + 1);
    const space_per_array_list = FOLD_OVERALLOCATE * (std.math.divCeil(u16, n_genes, n_arrays) catch unreachable);
    return struct {
        arraylists: [n_arrays]std.ArrayListUnmanaged(LastInt),

        const Self = @This();
        fn init(alloc: std.mem.Allocator) std.mem.Allocator.Error!@This() {
            var self = Self{ .arraylists = undefined };
            for (&self.arraylists) |*al| {
                al.* = try std.ArrayListUnmanaged(LastInt).initCapacity(alloc, space_per_array_list);
            }
            return self;
        }
        fn deinit(self: *Self, alloc: std.mem.Allocator) void {
            for (&self.arraylists) |*al| {
                al.deinit(alloc);
            }
        }
        // trys to append to set returns True if object already exists
        fn append(self: *Self, alloc: std.mem.Allocator, item: ByteArray(n_genes)) std.mem.Allocator.Error!bool {
            // check if item exists
            const first_byte = item.getByte();
            const last_int = item.getLastInt();
            const alp = &self.arraylists[first_byte];
            for (alp.items) |existing_int| {
                if (last_int == existing_int) {
                    return true;
                }
            }
            // insert in the item
            try alp.append(alloc, last_int);
            return false;
        }
        fn clear(self: *Self) void {
            for (&self.arraylists) |*al| {
                al.items.len = 0;
            }
        }
    };
}

test "init combo set" {
    var alloc = std.testing.allocator;
    var cs = try ComboSet(128).init(alloc);
    defer cs.deinit(alloc);
}

test "append combo set" {
    var alloc = std.testing.allocator;
    var cs = try ComboSet(100).init(alloc);
    defer cs.deinit(alloc);
    var bool_array = [_]bool{false} ++ [_]bool{true} ** 99;

    const ba = ByteArray(100).create(bool_array);
    var there = try cs.append(alloc, ba);
    try std.testing.expect(!there);
    there = try cs.append(alloc, ba);
    try std.testing.expect(there);
    try std.testing.expect(cs.arraylists[ba.getByte()].items.len == 1);
    bool_array = [_]bool{false} ++ [_]bool{true} ** 98 ++ [_]bool{false};
    there = try cs.append(alloc, ByteArray(100).create(bool_array));
    try std.testing.expect(!there);
    try std.testing.expect(cs.arraylists[ba.getByte()].items.len == 2);
}

test "clear combo set" {
    var alloc = std.testing.allocator;
    var cs = try ComboSet(100).init(alloc);
    defer cs.deinit(alloc);
    var bool_array = [_]bool{false} ++ [_]bool{true} ** 99;
    const ba = ByteArray(100).create(bool_array);
    var there = try cs.append(alloc, ba);
    try std.testing.expect(!there);
    there = try cs.append(alloc, ba);
    try std.testing.expect(there);
    cs.clear();
    there = try cs.append(alloc, ba);
    try std.testing.expect(!there);
}

test "clear array list" {
    var alloc = std.testing.allocator;
    var al = try std.ArrayListUnmanaged(u8).initCapacity(alloc, 10);
    defer al.deinit(alloc);
    try al.append(alloc, 10);
    print(al.items);
    al.items.len = 0;
    print(al.items);
    try al.append(alloc, 20);
    print(al.items);
}

/// return the number of collisions in one trial
fn oneTrial(comptime n_genes: u16, comptime n_types: usize, alloc: std.mem.Allocator, prob: f64, cs: *ComboSet(n_genes), prng: *std.rand.DefaultPrng) !u16 {
    // do random numbers
    const FloatVector = @Vector(n_types * n_genes, f64);
    const rand = prng.random();
    var rand_nums = try alloc.alloc(f64, n_types * n_genes);
    defer alloc.free(rand_nums);
    for (rand_nums) |*rand_num| {
        rand_num.* = rand.float(f64);
    }
    const bool_array: [n_types * n_genes]bool = @as(FloatVector, @splat(prob)) > @as(FloatVector, rand_nums[0 .. n_types * n_genes].*);
    // create a byte array for each gene being selected
    cs.clear();
    var first_index: usize = 0;
    var n_intersections: u16 = 0;
    while (first_index < rand_nums.len) : (first_index += n_genes) {
        const chunk: [n_genes]bool = bool_array[first_index..][0..n_genes].*;
        const there = try cs.append(alloc, ByteArray(n_genes).create(chunk));
        if (there) {
            n_intersections += 1;
        }
    }
    return n_intersections;
}

var test_prng = std.rand.DefaultPrng.init(0);
test "one trial" {
    var alloc = std.testing.allocator;
    const n_genes = 20;
    var cs = try ComboSet(n_genes).init(alloc);
    defer cs.deinit(alloc);
    var n_intersections = try oneTrial(n_genes, 10, alloc, 0, &cs, &test_prng);
    try std.testing.expect(n_intersections == 9);
    cs.clear();
    n_intersections = try oneTrial(n_genes, 5, alloc, 0.5, &cs, &test_prng);
    print(n_intersections);
    try std.testing.expect(n_intersections < 3);
}
