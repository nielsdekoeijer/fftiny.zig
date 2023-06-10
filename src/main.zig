const std = @import("std");
const fftiny = @import("fftiny.zig");
pub fn main() !void {
    const T = f32;
    const n: usize = 1024;

    @setEvalBranchQuota(10000000);
    var fft = comptime fftiny.plan(T, n, fftiny.FFTMODE.FW).init();
    var i = [_]std.math.Complex(T){
        std.math.Complex(T).init(0.5, 0.5),
        std.math.Complex(T).init(-0.5, -0.5),
    } ** (n / 2);
    var o = [_]std.math.Complex(T){std.math.Complex(T).init(0.0, 0.0)} ** (n);

    const p = 1000000;
    std.debug.print("Trying {} ffts of size {}", .{ p, n });
    for (1..p) |_| {
        fft.transform(&i, &o);
    }

    for (0..n) |q| {
        std.debug.print("{}", .{o[q]});
    }
}
