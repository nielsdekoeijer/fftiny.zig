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

pub fn initTwiddles(
    comptime T: type,
    N: usize,
    noalias twd: *[]const Complex(T),
) void {
    const M = @divExact(N, 8);
    const K = -2.0 * std.math.pi / @as(T, @floatFromInt(N));
    for (0..M) |i| {
        twd[i] = std.math.complex.exp(Complex(T).init(0.0, K * @as(T, @floatFromInt(i))));
    }
}

inline fn cpBF2(
    comptime T: type,
    h: usize,
    idx0: usize,
    idx1: usize,
    noalias inp: []const Complex(T),
    noalias out: *[]Complex(T),
) void {
    out.*[h + 0] = Complex(T).add(inp[idx0], inp[idx1]);
    out.*[h + 1] = Complex(T).sub(inp[idx0], inp[idx1]);
}

inline fn cpBF4(
    comptime T: type,
    comptime S: FFTDirection,
    l: usize,
    h: usize,
    k: usize,
    twd: Complex(T),
    noalias out: *[]Complex(T),
) void {
    const a = out.*[h + k + ((0 * l) / 4)];
    const b = out.*[h + k + ((1 * l) / 4)];
    const c = out.*[h + k + ((2 * l) / 4)];
    const d = out.*[h + k + ((3 * l) / 4)];

    switch (S) {
        FFTDirection.FW => blk: {
            const wc = Complex(T).mul(twd, c);
            const wd = Complex(T).mul(Complex(T).conjugate(twd), d);
            const p = Complex(T).add(wc, wd);
            const q = Complex(T).mulbyi(Complex(T).sub(wc, wd));

            out.*[h + k + ((l * 0) / 4)] = Complex(T).add(a, p);
            out.*[h + k + ((l * 1) / 4)] = Complex(T).sub(b, q);
            out.*[h + k + ((l * 2) / 4)] = Complex(T).sub(a, p);
            out.*[h + k + ((l * 3) / 4)] = Complex(T).add(b, q);
            break :blk;
        },

        FFTDirection.BW => blk: {
            const wc = Complex(T).mul(Complex(T).conjugate(twd), c);
            const wd = Complex(T).mul(twd, d);
            const p = Complex(T).add(wc, wd);
            const q = Complex(T).mulbyi(Complex(T).sub(wc, wd));

            out.*[h + k + ((l * 0) / 4)] = Complex(T).add(a, p);
            out.*[h + k + ((l * 1) / 4)] = Complex(T).add(b, q);
            out.*[h + k + ((l * 2) / 4)] = Complex(T).sub(a, p);
            out.*[h + k + ((l * 3) / 4)] = Complex(T).sub(b, q);
            break :blk;
        },
    }
}

pub fn cpFFT(
    comptime T: type,
    comptime S: FFTDirection,
    N: usize,
    noalias twd: []const Complex(T),
    noalias inp: []const Complex(T),
    noalias out: *[]Complex(T),
) void {
    const log2_n: usize = 31 - @clz(N);
    const r: usize = 32 - log2_n;
    var p: usize = 0;
    var q: usize = 0;

    var h2: usize = 0;
    for (0..N) |h| {
        h2 = h + 2;
        const c = 30 - @clz(h ^ h2);

        const _i0 = (p - q) >> r;
        const _i1 = _i0 ^ (N >> 1);
        
        if (c & 1) {
            out[h + 0] = inp[_i0];
            out[h + 1] = inp[_i1];
            cpBF4(1, out + h - 2, 1);
        } else {
        }
    }
}
