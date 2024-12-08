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

inline fn cpFFT(
    comptime T: type,
    comptime S: FFTDirection,
    comptime N: usize,
    noalias inp: []const Complex(T),
    noalias out: *[]Complex(T),
) void {
    comptime var p = 0;
    comptime var q = 0;
    inline for (0..N / 2) |g| {
        const c = comptime @floor(@log2(@as(f32, @floatFromInt(g ^ (g + 1)))));
        const idx0 = comptime @mod(@as(isize, @intCast(p - q)) / 4, N);
        const idx1 = comptime @mod(idx0 + N / 2, N);

        if (comptime @mod(c, 2) == 1) {
            out.*[comptime 2 * g + 0] = inp[comptime idx0];
            out.*[comptime 2 * g + 1] = inp[comptime idx1];
        } else {
            cpBF2(T, 2 * g, idx0, idx1, inp, out);
        } 

        if (comptime c > 1 - @mod(c, 2)) {
            comptime var j = 1 - @mod(c, 2);
            inline while (j < c) {
                inline for (0..comptime std.math.pow(usize, 2, j)) |b| {
                    cpBF4(
                        T,
                        S,
                        comptime std.math.pow(usize, 2, j + 2),
                        comptime 2 * g + 2 - std.math.pow(usize, 2, j + 2),
                        b,
                        out,
                    );
                }

                j += 2;
            }
        }

        const k = comptime @log2(@as(f32, @floatFromInt(N))) - c - 1;
        q = comptime @mod(q, std.math.pow(usize, 2, k)) + ((p >> k) & 1) * std.math.pow(usize, 2, k);
        p = comptime @mod(p, std.math.pow(usize, 2, k)) + (1 - (p >> k) & 1) * std.math.pow(usize, 2, k + 1);
    }
}

pub fn fwFFTIterative(
    comptime T: type,
    comptime n: usize,
    noalias inp: []const Complex(T),
    noalias out: *[]Complex(T),
) void {
    @setEvalBranchQuota(n * n * n);
    cpFFT(T, FFTDirection.FW, n, inp, out);
}

test "fwFFTIterative with 8-point FFT" {
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

    fwFFTIterative(T, 8, inps, @constCast(&outs));

    for (0..8) |i| {
        try std.testing.expectApproxEqAbs(ref[i].re, out[i].re, 1e-5);
        try std.testing.expectApproxEqAbs(ref[i].im, out[i].im, 1e-5);
    }
}

pub fn bwFFTIterative(comptime T: type, comptime n: usize, inp: []const Complex(T), out: *[]Complex(T)) void {
    @setEvalBranchQuota(n * n);
    cpFFT(T, FFTDirection.BW, n, inp, out);
}

test "bwFFTIterative with 8-point FFT" {
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

    bwFFTIterative(T, 8, inps, @constCast(&outs));

    for (0..8) |i| {
        try std.testing.expectApproxEqAbs(ref[i].re, out[i].re / 8.0, 1e-5);
        try std.testing.expectApproxEqAbs(ref[i].im, out[i].im / 8.0, 1e-5);
    }
}
