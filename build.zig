const std = @import("std");
const Build = std.Build;

pub fn build(b: *std.Build) void {
    const is_root = b.pkg_hash.len == 0;

    // these options are only applicable for local dev.
    // as a dependency, our module should inherit options from the one depending on it.
    const maybe_target = if (is_root) b.standardTargetOptions(.{}) else null;
    const maybe_optimize = if (is_root) b.standardOptimizeOption(.{}) else null;

    const install_step = b.getInstallStep();
    const unit_test_step = b.step("unit-test", "Run unit tests.");

    {
        const test_step = b.step("test", "Run all tests.");
        test_step.dependOn(unit_test_step);
    }

    {
        const check_step = b.step("check", "Check step.");
        check_step.dependOn(install_step);
    }

    const objwf_mod = b.addModule("objwf", .{
        .root_source_file = b.path("src/objwf.zig"),
        .target = maybe_target,
        .optimize = maybe_optimize,
    });

    // the rest of this code is only applicable for local dev.
    if (!is_root) return;
    const bin_opts: BinOptions = .fromBuildOptions(b);
    const filters = b.option([]const []const u8, "filter", "List of filters for tests.") orelse &.{};

    const unit_test_exe = b.addTest(.{
        .name = "unit-test",
        .root_module = objwf_mod,
        .filters = filters,
    });
    const unit_test_output = addExeOutputs(b, .{
        .exe = unit_test_exe,
        .step = unit_test_step,
        .bin = bin_opts,
        .install = .{},
    });
    _ = unit_test_output;
}

const BinOptions = packed struct {
    install: bool,
    run: bool,

    fn fromBuildOptions(b: *Build) BinOptions {
        const no_bin = b.option(
            bool,
            "no-bin",
            "Don't install any of the binaries implied by the specified steps.",
        ) orelse false;
        const no_run = b.option(
            bool,
            "no-run",
            "Don't run any of the executables implied by the specified steps.",
        ) orelse false;
        return .{
            .install = !no_bin,
            .run = !no_run,
        };
    }
};

const ExeOutputs = struct {
    install: ?*Build.Step.InstallArtifact,
    run: ?*Build.Step.Run,
};

fn addExeOutputs(
    b: *Build,
    params: struct {
        exe: *Build.Step.Compile,
        step: *Build.Step,
        bin: BinOptions,
        install: Build.Step.InstallArtifact.Options,
    },
) ExeOutputs {
    const exe = params.exe;
    const step = params.step;
    const bin_opts = params.bin;
    const install_opts = params.install;

    const maybe_exe_install = if (bin_opts.install) b.addInstallArtifact(exe, install_opts) else null;
    const maybe_exe_run = if (bin_opts.run) b.addRunArtifact(exe) else null;

    step.dependOn(&exe.step);
    if (maybe_exe_install) |exe_install| step.dependOn(&exe_install.step);
    if (maybe_exe_run) |exe_run| step.dependOn(&exe_run.step);

    const install_step = b.getInstallStep();
    install_step.dependOn(&exe.step);
    if (maybe_exe_install) |exe_install| install_step.dependOn(&exe_install.step);

    return .{
        .install = maybe_exe_install,
        .run = maybe_exe_run,
    };
}
