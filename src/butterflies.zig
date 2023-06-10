const std = @import("std");
const utility = @import("utility.zig");

pub const FFTMODE = enum { FW, BW };

// =butterfly-2================================================================
// strided list
pub inline fn butterfly2sl(comptime T: type, comptime mode: FFTMODE, x0: [2]T, x1: [2]T) [2][2]T {
    _ = mode;
    const x0v = @Vector(2, T){ x0[0], x0[1] };
    const x1v = @Vector(2, T){ x1[0], x1[1] };
    return [_][2]T{ x0v + x1v, x0v - x1v };
}

test "butterfly2sl" {
    const T = f32;
    var x0: @Vector(2, T) = [_]T{
        -2.3741363830674853,
        -0.27588287433592273,
    };
    var x1: @Vector(2, T) = [_]T{
        -2.3846697075321135,
        -2.6767359305510112,
    };
    const y = butterfly2sl(T, FFTMODE.FW, x0, x1);
    const r: [2]@Vector(2, T) = [_]@Vector(2, T){ @Vector(2, T){
        -4.758806090599599,
        -2.952618804886934,
    }, @Vector(2, T){
        0.01053332446462818,
        2.4008530562150883,
    } };

    for (0..2) |i| {
        try std.testing.expect(std.math.approxEqAbs(T, y[0][i], r[0][i], 1e-6));
        try std.testing.expect(std.math.approxEqAbs(T, y[1][i], r[1][i], 1e-6));
    }
}

// strided packed
pub inline fn butterfly2sp(comptime T: type, comptime mode: FFTMODE, x0: [2]T, x1: [2]T) [4]T {
    _ = mode;
    return [_]T{
        x0[0] + x1[0],
        x0[1] + x1[1],
        x0[0] - x1[0],
        x0[1] - x1[1],
    };
}

test "butterfly2sp" {
    const T = f32;
    var x0: @Vector(2, T) = [_]T{
        -2.3741363830674853,
        -0.27588287433592273,
    };

    var x1: @Vector(2, T) = [_]T{
        -2.3846697075321135,
        -2.6767359305510112,
    };

    const y = butterfly2sp(T, FFTMODE.FW, x0, x1);

    const r: @Vector(4, T) = [_]T{
        -4.758806090599599,
        -2.952618804886934,
        0.01053332446462818,
        2.4008530562150883,
    };

    for (0..4) |i| {
        try std.testing.expect(std.math.approxEqAbs(T, y[i], r[i], 1e-6));
    }
}

// 2x interleaved list
pub inline fn butterfly2x2il(comptime T: type, comptime mode: FFTMODE, x02: [4]T, x13: [4]T) [2][4]T {
    _ = mode;
    const x02v = @Vector(4, T){ x02[0], x02[1], x02[2], x02[3] };
    const x13v = @Vector(4, T){ x13[0], x13[1], x13[2], x13[3] };
    return [_][4]T{ x02v + x13v, x02v - x13v };
}

test "butterfly2x2il" {
    const T = f32;
    var x0: @Vector(4, T) = [_]T{
        1.0,
        2.0,
        3.0,
        4.0,
    };

    var x1: @Vector(4, T) = [_]T{
        5.0,
        6.0,
        7.0,
        8.0,
    };

    const y = butterfly2x2il(T, FFTMODE.FW, x0, x1);

    const r0: @Vector(4, T) = [_]T{
        1.0 + 5.0,
        2.0 + 6.0,
        3.0 + 7.0,
        4.0 + 8.0,
    };

    const r1: @Vector(4, T) = [_]T{
        1.0 - 5.0,
        2.0 - 6.0,
        3.0 - 7.0,
        4.0 - 8.0,
    };

    for (0..4) |i| {
        try std.testing.expect(std.math.approxEqAbs(T, y[0][i], r0[i], 1e-6));
        try std.testing.expect(std.math.approxEqAbs(T, y[1][i], r1[i], 1e-6));
    }
}

// 2x interleaved packed
pub inline fn butterfly2x2ip(comptime T: type, comptime mode: FFTMODE, x02: [4]T, x13: [4]T) [8]T {
    _ = mode;
    return [8]T{
        x02[0] + x13[0],
        x02[1] + x13[1],
        x02[2] + x13[2],
        x02[3] + x13[3],
        x02[0] - x13[0],
        x02[1] - x13[1],
        x02[2] - x13[2],
        x02[3] - x13[3],
    };
}

test "butterfly2x2ip" {
    const T = f32;
    var x0: @Vector(4, T) = [_]T{
        1.0,
        2.0,
        3.0,
        4.0,
    };

    var x1: @Vector(4, T) = [_]T{
        5.0,
        6.0,
        7.0,
        8.0,
    };

    const y = butterfly2x2ip(T, FFTMODE.FW, x0, x1);

    const r: @Vector(8, T) = [_]T{
        1.0 + 5.0,
        2.0 + 6.0,
        3.0 + 7.0,
        4.0 + 8.0,
        1.0 - 5.0,
        2.0 - 6.0,
        3.0 - 7.0,
        4.0 - 8.0,
    };

    for (0..8) |i| {
        try std.testing.expect(std.math.approxEqAbs(T, y[i], r[i], 1e-6));
    }
}

// 2x2 strided list
pub inline fn butterfly2x2sl(comptime T: type, comptime mode: FFTMODE, x01: [4]T, x23: [4]T) [2][4]T {
    const temp = utility.transpose2x2(T, mode, x01, x23);
    return butterfly2x2il(T, mode, temp[0], temp[1]);
}

test "butterfly2x2sl" {
    const T = f32;
    var x0: @Vector(4, T) = [_]T{
        1.0,
        2.0,
        5.0,
        6.0,
    };

    var x1: @Vector(4, T) = [_]T{
        3.0,
        4.0,
        7.0,
        8.0,
    };

    const y = butterfly2x2sl(T, FFTMODE.FW, x0, x1);

    const r0: @Vector(4, T) = [_]T{
        1.0 + 5.0,
        2.0 + 6.0,
        3.0 + 7.0,
        4.0 + 8.0,
    };

    const r1: @Vector(4, T) = [_]T{
        1.0 - 5.0,
        2.0 - 6.0,
        3.0 - 7.0,
        4.0 - 8.0,
    };

    for (0..4) |i| {
        try std.testing.expect(std.math.approxEqAbs(T, y[0][i], r0[i], 1e-6));
        try std.testing.expect(std.math.approxEqAbs(T, y[1][i], r1[i], 1e-6));
    }
}

// 2x2 strided packed
pub inline fn butterfly2x2sp(comptime T: type, comptime mode: FFTMODE, x01: [4]T, x23: [4]T) [8]T {
    const temp = utility.transpose2x2(T, mode, x01, x23);
    return butterfly2x2ip(T, mode, temp[0], temp[1]);
}

test "butterfly2x2sp" {
    const T = f32;
    var x0: @Vector(4, T) = [_]T{
        1.0,
        2.0,
        5.0,
        6.0,
    };

    var x1: @Vector(4, T) = [_]T{
        3.0,
        4.0,
        7.0,
        8.0,
    };

    const y = butterfly2x2sp(T, FFTMODE.FW, x0, x1);

    const r: @Vector(8, T) = [_]T{
        1.0 + 5.0,
        2.0 + 6.0,
        3.0 + 7.0,
        4.0 + 8.0,
        1.0 - 5.0,
        2.0 - 6.0,
        3.0 - 7.0,
        4.0 - 8.0,
    };

    for (0..8) |i| {
        try std.testing.expect(std.math.approxEqAbs(T, y[i], r[i], 1e-6));
    }
}

// =butterfly-4================================================================
// strided packed
pub inline fn butterfly4sp(comptime T: type, comptime mode: FFTMODE, x01: [4]T, x23: [4]T) [8]T {
    var y = butterfly2x2il(T, mode, x01, x23);

    y[1][2..4].* = utility.rotateN90(T, mode, y[1][2..4].*);

    return butterfly2x2sp(T, mode, y[0], y[1]);
}

test "butterfly4sp" {
    const T = f32;

    var x0: @Vector(4, T) = [_]T{
        1.2514228558126792,
        0.2574792970535659,
        -1.002327101854918,
        -1.1751269599415906,
    };

    var x1: @Vector(4, T) = [_]T{
        -1.3050896806819041,
        -1.30529768383336,
        0.22785218401785495,
        0.11806180009255686,
    };

    const y = butterfly4sp(T, FFTMODE.FW, x0, x1);

    const r: @Vector(8, T) = [_]T{
        -0.828141742706288,
        -2.1048835466288276,
        1.2633237764604361,
        2.792956266759699,
        0.7208080929678382,
        0.0092467730692396,
        3.849701296528731,
        0.33259769501415315,
    };

    for (0..8) |i| {
        try std.testing.expect(std.math.approxEqAbs(T, y[i], r[i], 1e-6));
    }
}

// =butterfly-sr-4=============================================================
pub inline fn butterflySR4l(comptime T: type, comptime mode: FFTMODE, comptime size: usize, comptime k: usize, x0: *[2]T, x1: *[2]T, x2: *[2]T, x3: *[2]T) void {
    const twr = comptime utility.twiddle(T, mode, size, k).re;
    const twi = comptime utility.twiddle(T, mode, size, k).im;

    const tw1 = comptime @Vector(2, T){ twr, twi };
    const tw3 = comptime @Vector(2, T){ twr, -twi };

    // TODO: parallize this ?
    const odds1 = cmul(T, tw1, x1.*);
    const odds3 = cmul(T, tw3, x3.*);

    var odds = butterfly2sl(T, mode, odds1, odds3);
    odds[1] = utility.rotateN90(T, mode, odds[1]);

    const a = butterfly2sl(T, mode, x0.*, odds[0]);
    const b = butterfly2sl(T, mode, x2.*, odds[1]);

    x0.* = a[0];
    x1.* = a[1];
    x2.* = b[0];
    x3.* = b[1];
}
