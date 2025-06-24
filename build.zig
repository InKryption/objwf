const std = @import("std");
const Build = std.Build;

pub fn build(b: *std.Build) void {
    const is_root = b.pkg_hash.len == 0;
    const target = if (is_root) b.standardTargetOptions(.{}) else null;
    const optimize = if (is_root) b.standardOptimizeOption(.{}) else null;

    const bin_install = !(b.option(bool, "no-bin", "Don't install any of the binaries implied by the specified steps.") orelse false);
    const bin_run = !(b.option(bool, "no-run", "Don't run any of the executables implied by the specified steps.") orelse false);
    const filters = b.option([]const []const u8, "filter", "List of filters for tests.") orelse &.{};

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
        .target = target,
        .imports = &.{},
    });

    if (is_root) {
        const unit_test_exe = b.addTest(.{
            .name = "unit-test",
            .root_module = objwf_mod,
            .optimize = optimize orelse .Debug,
            .filters = filters,
        });
        const unit_test_output = addInstallAndRun(b, unit_test_exe, unit_test_step, .{
            .run = bin_run,
            .install = bin_install,
            .install_options = .{},
        });
        _ = unit_test_output;
    }
}

const InstallAndRunOptions = struct {
    run: bool,
    install: bool,
    install_options: Build.Step.InstallArtifact.Options,
};
const InstallAndRun = struct {
    install: ?*Build.Step.InstallArtifact,
    run: ?*Build.Step.Run,
};

fn addInstallAndRun(
    b: *Build,
    exe: *Build.Step.Compile,
    step: *Build.Step,
    options: InstallAndRunOptions,
) InstallAndRun {
    const install_step = b.getInstallStep();
    const maybe_exe_install = if (options.install) b.addInstallArtifact(exe, options.install_options) else null;
    const maybe_exe_run = if (options.run) b.addRunArtifact(exe) else null;

    step.dependOn(&exe.step);
    if (maybe_exe_install) |exe_install| step.dependOn(&exe_install.step);
    if (maybe_exe_run) |exe_run| step.dependOn(&exe_run.step);

    install_step.dependOn(&exe.step);
    if (maybe_exe_install) |exe_install| install_step.dependOn(&exe_install.step);

    return .{
        .install = maybe_exe_install,
        .run = maybe_exe_run,
    };
}
