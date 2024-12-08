const std = @import("std");
const testing = std.testing;
const Complex = std.math.complex.Complex;

const FFTDirection = enum {
    FW,
    BW,
};

const FFTNormalization = enum {
    FW,
    BW,
    ORTHO,
};

inline fn cpBF2(
    comptime T: type,
    comptime h: usize,
    comptime idx0: usize,
    comptime idx1: usize,
    noalias inp: []const Complex(T),
    noalias out: *[]Complex(T),
) void {
    out.*[h + 0] = Complex(T).add(inp[idx0], inp[idx1]);
    out.*[h + 1] = Complex(T).sub(inp[idx0], inp[idx1]);
}

inline fn cpBF4(
    comptime T: type,
    comptime S: FFTDirection,
    comptime l: usize,
    comptime h: usize,
    comptime k: usize,
    noalias out: *[]Complex(T),
) void {
    const a = out.*[comptime h + k + ((0 * l) / 4)];
    const b = out.*[comptime h + k + ((1 * l) / 4)];
    const c = out.*[comptime h + k + ((2 * l) / 4)];
    const d = out.*[comptime h + k + ((3 * l) / 4)];
    const w = comptime std.math.complex.exp(Complex(T).init(0.0, -2.0 * std.math.pi * @as(T, @floatFromInt(k)) / @as(T, @floatFromInt(l))));
    const m = comptime Complex(T).conjugate(w);

    switch (S) {
        FFTDirection.FW => blk: {
            // const wc = Complex(T).mul(w, c);
            // const wd = Complex(T).mul(comptime Complex(T).conjugate(w), d);
            const wm = @Vector(4, T) { w.re, w.im, m.re, m.im };
            const cd = @Vector(4, T) { c.re, c.im, d.re, d.im };
            const re = @shuffle(f32, wm, wm, @Vector(4, i32){ 0, 0, 2, 2 }) * cd;
            const _b = @shuffle(f32, cd, cd, @Vector(4, i32){ 1, 0, 3, 2 });
            const im = @shuffle(f32, wm, wm, @Vector(4, i32){ 1, 1, 3, 3 }) * _b;
            const wcmd = re + @Vector(4, f32){ -im[0], im[1], -im[2], im[3] };

            // const p = Complex(T).add(wc, wd);
            // const q = Complex(T).mulbyi(Complex(T).sub(wc, wd));
            const wcwc = @shuffle(f32, wcmd, undefined, @Vector(4, i32) {0, 1, 0, 1});
            const mdmd = @shuffle(f32, wcmd, -wcmd, @Vector(4, i32) {2, 3, -3, -4});
            const ab = @Vector(4, T) { a.re, a.im, b.re, b.im };
            var pq = wcwc + mdmd;
            pq = @Vector(4, f32) {pq[0], pq[1], pq[3], -pq[2] };

            {
                const rs = ab + pq;
                out.*[comptime h + k + ((l * 0) / 4)] = Complex(T).init(rs[0], rs[1]);
                out.*[comptime h + k + ((l * 1) / 4)] = Complex(T).init(rs[2], rs[3]);
            }
            {
                const rs = ab - pq;
                out.*[comptime h + k + ((l * 2) / 4)] = Complex(T).init(rs[0], rs[1]);
                out.*[comptime h + k + ((l * 3) / 4)] = Complex(T).init(rs[2], rs[3]);
            }

            break :blk;
        },

        FFTDirection.BW => blk: {
            const wc = Complex(T).mul(comptime Complex(T).conjugate(w), c);
            const wd = Complex(T).mul(w, d);
            const p = Complex(T).add(wc, wd);
            const q = Complex(T).mulbyi(Complex(T).sub(wc, wd));

            out.*[comptime h + k + ((l * 0) / 4)] = Complex(T).add(a, p);
            out.*[comptime h + k + ((l * 1) / 4)] = Complex(T).add(b, q);
            out.*[comptime h + k + ((l * 2) / 4)] = Complex(T).sub(a, p);
            out.*[comptime h + k + ((l * 3) / 4)] = Complex(T).sub(b, q);
            break :blk;
        },
    }
}

pub fn cpFFTRecursive(
    comptime T: type,
    comptime S: FFTDirection,
    comptime n: usize,
    comptime l: usize,
    comptime h: usize,
    comptime i: isize,
    noalias inp: []const Complex(T),
    noalias out: *[]Complex(T),
) void {
    const idx0 = @as(usize, @intCast(@mod(i, @as(isize, @intCast(n)))));
    switch (l) {
        1 => blk: {
            out.*[h] = inp[idx0];
            break :blk;
        },
        2 => blk: {
            const idx1 = @as(usize, @intCast(@mod(i + @as(isize, @intCast(n / l)), @as(isize, @intCast(n)))));
            cpBF2(T, h, idx0, idx1, inp, out);
            break :blk;
        },
        else => blk: {
            cpFFTRecursive(T, S, n, l / 2, h, i, inp, out);
            cpFFTRecursive(T, S, n, l / 4, h + ((1 * l) / 2), i + @as(isize, @intCast(n / l)), inp, out);
            cpFFTRecursive(T, S, n, l / 4, h + ((3 * l) / 4), i - @as(isize, @intCast(n / l)), inp, out);
            inline for (0..l / 4) |b| {
                cpBF4(T, S, l, h, b, out);
            }
            break :blk;
        },
    }
}

pub fn fwFFTRecursive(comptime T: type, comptime n: usize, inp: []const Complex(T), out: *[]Complex(T)) void {
    @setEvalBranchQuota(n * n);
    cpFFTRecursive(T, FFTDirection.FW, n, n, 0, 0, inp, out);
}

test "fwFFT with 8-point FFT" {
    const T = f32;

    var inp = [8]Complex(T){
        Complex(T){ .re = -0.71134056, .im = -0.34246695 },
        Complex(T){ .re = -0.34408238, .im = 0.20285457 },
        Complex(T){ .re = -0.72280309, .im = 0.79882934 },
        Complex(T){ .re = -1.38351121, .im = 0.13329499 },
        Complex(T){ .re = -1.07962283, .im = 0.57833735 },
        Complex(T){ .re = 1.85331861, .im = -0.00671944 },
        Complex(T){ .re = 0.10836201, .im = -0.10296904 },
        Complex(T){ .re = 0.41946510, .im = 0.03223966 },
    };

    const ref = [8]Complex(T){
        Complex(T){ .re = -1.86021434, .im = 1.29340048 },
        Complex(T){ .re = 1.21082841, .im = 2.81578900 },
        Complex(T){ .re = -1.14592182, .im = -2.93327223 },
        Complex(T){ .re = -0.03496764, .im = 0.99999023 },
        Complex(T){ .re = -2.95059458, .im = 0.57006091 },
        Complex(T){ .re = 1.32933291, .im = -2.99506740 },
        Complex(T){ .re = -1.20712279, .im = 2.01329244 },
        Complex(T){ .re = -1.03206459, .im = -4.50392900 },
    };

    var out: [8]Complex(T) = [8]Complex(T){
        Complex(T){ .re = 0.0, .im = 0.0 },
        Complex(T){ .re = 0.0, .im = 0.0 },
        Complex(T){ .re = 0.0, .im = 0.0 },
        Complex(T){ .re = 0.0, .im = 0.0 },
        Complex(T){ .re = 0.0, .im = 0.0 },
        Complex(T){ .re = 0.0, .im = 0.0 },
        Complex(T){ .re = 0.0, .im = 0.0 },
        Complex(T){ .re = 0.0, .im = 0.0 },
    };

    const inps: []const Complex(T) = inp[0..];
    const outs: []Complex(T) = out[0..];

    fwFFTRecursive(T, 8, inps, @constCast(&outs));

    for (0..8) |i| {
        try std.testing.expectApproxEqAbs(ref[i].re, out[i].re, 1e-5);
        try std.testing.expectApproxEqAbs(ref[i].im, out[i].im, 1e-5);
    }
}

pub fn bwFFTRecursive(comptime T: type, comptime n: usize, inp: []const Complex(T), out: *[]Complex(T)) void {
    @setEvalBranchQuota(n * n);
    cpFFTRecursive(T, FFTDirection.BW, n, n, 0, 0, inp, out);
}

test "bwFFT with 8-point FFT" {
    const T = f32;

    var inp = [8]Complex(T){
        Complex(T){ .re = -1.86021434, .im = 1.29340048 },
        Complex(T){ .re = 1.21082841, .im = 2.81578900 },
        Complex(T){ .re = -1.14592182, .im = -2.93327223 },
        Complex(T){ .re = -0.03496764, .im = 0.99999023 },
        Complex(T){ .re = -2.95059458, .im = 0.57006091 },
        Complex(T){ .re = 1.32933291, .im = -2.99506740 },
        Complex(T){ .re = -1.20712279, .im = 2.01329244 },
        Complex(T){ .re = -1.03206459, .im = -4.50392900 },
    };

    const ref = [8]Complex(T){
        Complex(T){ .re = -0.71134056, .im = -0.34246695 },
        Complex(T){ .re = -0.34408238, .im = 0.20285457 },
        Complex(T){ .re = -0.72280309, .im = 0.79882934 },
        Complex(T){ .re = -1.38351121, .im = 0.13329499 },
        Complex(T){ .re = -1.07962283, .im = 0.57833735 },
        Complex(T){ .re = 1.85331861, .im = -0.00671944 },
        Complex(T){ .re = 0.10836201, .im = -0.10296904 },
        Complex(T){ .re = 0.41946510, .im = 0.03223966 },
    };

    var out: [8]Complex(T) = [8]Complex(T){
        Complex(T){ .re = 0.0, .im = 0.0 },
        Complex(T){ .re = 0.0, .im = 0.0 },
        Complex(T){ .re = 0.0, .im = 0.0 },
        Complex(T){ .re = 0.0, .im = 0.0 },
        Complex(T){ .re = 0.0, .im = 0.0 },
        Complex(T){ .re = 0.0, .im = 0.0 },
        Complex(T){ .re = 0.0, .im = 0.0 },
        Complex(T){ .re = 0.0, .im = 0.0 },
    };

    const inps: []const Complex(T) = inp[0..];
    const outs: []Complex(T) = out[0..];

    bwFFTRecursive(T, 8, inps, @constCast(&outs));

    for (0..8) |i| {
        try std.testing.expectApproxEqAbs(ref[i].re, out[i].re / 8.0, 1e-5);
        try std.testing.expectApproxEqAbs(ref[i].im, out[i].im / 8.0, 1e-5);
    }
}
