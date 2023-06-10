// =utils======================================================================
// transpose 2x2
pub inline fn transpose2x2(comptime T: type, comptime mode: FFTMODE, x01: [4]T, x23: [4]T) [2][4]T {
    _ = mode;
    return [_][4]T{
        [_]T{ x01[0], x01[1], x23[0], x23[1] },
        [_]T{ x01[2], x01[3], x23[2], x23[3] },
    };
}

test "transpose 2x2" {
    const T = f32;
    const x01: @Vector(4, T) = [_]T{ 1.0, 2.0, 3.0, 4.0 };
    const x23: @Vector(4, T) = [_]T{ 5.0, 6.0, 7.0, 8.0 };

    const y = transpose2x2(T, FFTMODE.FW, x01, x23);

    const x02: @Vector(4, T) = [_]T{ 1.0, 2.0, 5.0, 6.0 };
    const x13: @Vector(4, T) = [_]T{ 3.0, 4.0, 7.0, 8.0 };
    for (0..4) |i| {
        try std.testing.expect(y[0][i] == x02[i]);
        try std.testing.expect(y[1][i] == x13[i]);
    }
}

// rotate 90 aka mul by -i
pub inline fn rotateN90(comptime T: type, comptime mode: FFTMODE, x: [2]T) [2]T {
    switch (comptime mode) {
        FFTMODE.FW => return [_]T{ x[1], -x[0] },
        FFTMODE.BW => return [_]T{ -x[1], x[0] },
    }
}

test "rotate 90" {
    const T = f32;
    const x01: @Vector(2, T) = [_]T{ 1.0, 2.0 };
    const y01 = rotateN90(T, FFTMODE.FW, x01);
    const z01: @Vector(2, T) = [_]T{ 2.0, -1.0 };
    try std.testing.expect(y01[0] == z01[0]);
    try std.testing.expect(y01[1] == z01[1]);
}

// complex mul
pub inline fn cmul(comptime T: type, comptime x: [2]T, y: [2]T) [2]T {
    if (comptime (x[1] == 0.0)) {
        // ergo x[0] == + / 1 1.0...
        if (x[0] > 0.0) {
            return y;
        } else {
            return [_]T{ -y[0], -y[1] };
        }
    }

    if (comptime (x[0] == 0.0)) {
        // ergo x[1] == + / 1 1.0...
        if (x[1] > 0.0) {
            return [_]T{ -y[1], y[0] };
        } else {
            return [_]T{ y[1], -y[0] };
        }
        return [_]T{ -x[1] * y[1], x[1] * y[0] };
    }

    const temp = @Vector(4, T){ x[0], x[0], x[1], -x[1] } * @Vector(4, T){ y[0], y[1], y[0], y[1] };
    return [_]T{ temp[0] + temp[3], temp[1] + temp[2] };
}

test "cmul" {
    const T = f32;
    const x: @Vector(2, T) = comptime [_]T{ 1.0, 2.0 };
    const y: @Vector(2, T) = [_]T{ 3.0, 4.0 };
    const z: @Vector(2, T) = cmul(T, x, y);
    const r: @Vector(2, T) = [_]T{ -5.0, 10.0 };
    try std.testing.expect(z[0] == r[0]);
    try std.testing.expect(z[1] == r[1]);
}

pub inline fn twiddle(comptime T: type, comptime mode: FFTMODE, comptime size: usize, comptime k: usize) std.math.Complex(T) {
    const phase = switch (comptime mode) {
        FFTMODE.FW => -2.0 * std.math.pi / (4 * @as(T, size)) * @as(T, k),
        FFTMODE.BW => 2.0 * std.math.pi / (4 * @as(T, size)) * @as(T, k),
    };

    const cos = @cos(phase);
    const sin = @sin(phase);

    return std.math.Complex(T).init(cos, sin);
}

test "twiddle" {}

