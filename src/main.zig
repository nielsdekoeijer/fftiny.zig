const std = @import("std");
const Complex = std.math.complex.Complex;

pub fn main() !void {
    const T = f32;

    const fwFFTIterative = @import("iterative.zig").fwFFTIterative;

    var inp = try std.heap.page_allocator.alloc(Complex(T), 16);
    defer std.heap.page_allocator.free(inp);

    var out = try std.heap.page_allocator.alloc(Complex(T), 16);
    defer std.heap.page_allocator.free(out);

    @breakpoint();
    fwFFTIterative(T, comptime 16, inp[0..], @constCast(&out[0..]));
}
