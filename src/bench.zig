const std = @import("std");
const Complex = std.math.complex.Complex;

pub fn fwFFTIterativeBenchmark(comptime T: type, comptime niter: usize, comptime size: usize, stdout: anytype, bw: anytype) !void {
    const fwFFTIterative = @import("iterative.zig").fwFFTIterative;

    var inp = try std.heap.page_allocator.alloc(Complex(T), size);
    defer std.heap.page_allocator.free(inp);

    var out = try std.heap.page_allocator.alloc(Complex(T), size);
    defer std.heap.page_allocator.free(out);

    const fw_start = std.time.nanoTimestamp();
    for (0..niter) |_| {
        fwFFTIterative(T, comptime size, inp[0..], @constCast(&out[0..]));
    }
    const fw_end = std.time.nanoTimestamp();

    try stdout.print(
        "-> {s:20}: {d:10.2} ns, mean: {d:10.2} ns\n",
        .{ "fwFFTIterative", fw_end - fw_start, @as(T, @floatFromInt(fw_end - fw_start)) / @as(T, @floatFromInt(niter)) },
    );
    try bw.flush();
}

pub fn fwFFTRecursiveBenchmark(comptime T: type, comptime niter: usize, comptime size: usize, stdout: anytype, bw: anytype) !void {
    const fwFFTRecursive = @import("recursive.zig").fwFFTRecursive;

    var inp = try std.heap.page_allocator.alloc(Complex(T), size);
    defer std.heap.page_allocator.free(inp);

    var out = try std.heap.page_allocator.alloc(Complex(T), size);
    defer std.heap.page_allocator.free(out);

    const fw_start = std.time.nanoTimestamp();
    for (0..niter) |_| {
        fwFFTRecursive(T, comptime size, inp[0..], @constCast(&out[0..]));
    }
    const fw_end = std.time.nanoTimestamp();

    try stdout.print(
        "-> {s:20}: {d:10.2} ns, mean: {d:10.2} ns\n",
        .{ "fwFFTRecursive", fw_end - fw_start, @as(T, @floatFromInt(fw_end - fw_start)) / @as(T, @floatFromInt(niter)) },
    );
    try bw.flush();
}

pub fn fwFFTWBenchmark(comptime T: type, comptime niter: usize, comptime size: usize, stdout: anytype, bw: anytype) !void {
    const fftw = @cImport(@cInclude("fftw3.h"));
    const fftwf_complex = [2]f32;

    const inp_ptr = fftw.fftwf_malloc(@sizeOf(fftwf_complex) * size) orelse return error.OutOfMemory;
    defer fftw.fftwf_free(inp_ptr);
    const inp = @as([*c][2]f32, @alignCast(@ptrCast(inp_ptr)));

    const out_ptr = fftw.fftwf_malloc(@sizeOf(fftwf_complex) * size) orelse return error.OutOfMemory;
    defer fftw.fftwf_free(out_ptr);
    const out = @as([*c][2]f32, @alignCast(@ptrCast(out_ptr)));

    const plan = fftw.fftwf_plan_dft_1d(@intCast(size), inp, out, fftw.FFTW_FORWARD, fftw.FFTW_ESTIMATE);
    if (plan == null) return error.PlanCreationFailed;
    defer fftw.fftwf_destroy_plan(plan);

    const fw_start = std.time.nanoTimestamp();
    for (0..niter) |_| {
        fftw.fftwf_execute(plan);
    }
    const fw_end = std.time.nanoTimestamp();

    // Output results to stdout
    try stdout.print(
        "-> {s:20}: {d:10.2} ns, mean: {d:10.2} ns\n",
        .{ "fwFFTW", fw_end - fw_start, @as(T, @floatFromInt(fw_end - fw_start)) / @as(T, @floatFromInt(niter)) },
    );
    try bw.flush();
}

pub fn main() !void {
    const sizes = [_]usize{ 8, 16, 32, 64, 128 };
    const niter = 10000000;
    const T = f32;

    const stdout_file = std.io.getStdOut().writer();
    var bw = std.io.bufferedWriter(stdout_file);
    const stdout = bw.writer();

    inline for (sizes) |size| {
        try stdout.print("Starting size {}...\n", .{size});
        try bw.flush();

        try fwFFTIterativeBenchmark(T, niter, comptime size, &stdout, &bw);
        try fwFFTRecursiveBenchmark(T, niter, comptime size, &stdout, &bw);
        try fwFFTWBenchmark(T, niter, comptime size, &stdout, &bw);

        try stdout.print("\n", .{});
        try bw.flush();
    }

    try bw.flush();
}
