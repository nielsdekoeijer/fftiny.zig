const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const lib = b.addStaticLibrary(.{
        .name = "fftiny",
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });
    b.installArtifact(lib);

    const bench_exe = b.addExecutable(.{
        .name = "fftiny",
        .root_source_file = b.path("src/bench.zig"),
        .target = target,
        .optimize = std.builtin.OptimizeMode.ReleaseFast,
    });
    bench_exe.linkSystemLibrary("fftw3f"); // for single precision
    bench_exe.linkLibC();
    b.installArtifact(bench_exe);

    const bench_run_cmd = b.addRunArtifact(bench_exe);
    bench_run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        bench_run_cmd.addArgs(args);
    }
    const bench_run_step = b.step("bench", "Do the bench");
    bench_run_step.dependOn(&bench_run_cmd.step);

    const lib_unit_tests = b.addTest(.{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });
    const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);
}
