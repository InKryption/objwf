const std = @import("std");
const stridx = @import("stridx.zig");

comptime {
    _ = stridx;
}

const StrSet = @This();
bytes: []const u8,
indexer: Indexer,

pub const Index = stridx.Index;
pub const Indexer = stridx.Indexer;

pub const empty: StrSet = .{
    .bytes = &.{},
    .set = .empty,
};

pub fn deinit(self: StrSet, gpa: std.mem.Allocator) void {
    gpa.free(self.bytes);
    var set = self.indexer;
    set.deinit(gpa);
}

pub fn getStr(self: StrSet, index: StrSet.Index) []const u8 {
    std.debug.assert(self.indexer.containsContext(index, stridx.hashCtx(self.bytes)));
    return stridx.idxSliceTo(self.bytes, index);
}

pub fn getIndex(self: StrSet, str: []const u8) ?StrSet.Index {
    return self.indexer.getKeyAdapted(str, stridx.hashStrAdapter(self.bytes));
}

pub fn fmt(self: StrSet, index: StrSet.Index) Fmt {
    return .{
        .index = index,
        .str_set = self,
    };
}

pub const Fmt = struct {
    index: StrSet.Index,
    str_set: StrSet,

    pub fn format(
        self: Fmt,
        comptime fmt_str: []const u8,
        fmt_options: std.fmt.FormatOptions,
        writer: anytype,
    ) @TypeOf(writer).Error!void {
        try std.fmt.formatType(
            self.str_set.getStr(self.index),
            fmt_str,
            fmt_options,
            writer,
            std.options.fmt_max_depth,
        );
    }
};

pub const ListFmt = struct {
    indices: []const StrSet.Index,
    str_set: StrSet,

    pub fn format(
        self: ListFmt,
        comptime fmt_str: []const u8,
        fmt_options: std.fmt.FormatOptions,
        writer: anytype,
    ) @TypeOf(writer).Error!void {
        for (self.indices, 0..) |index, i| {
            if (i != 0) try writer.writeAll(" ");
            try index.fmt(self.str_set).format(fmt_str, fmt_options, writer);
        }
    }
};

pub fn listFmt(str_set: StrSet, indices: []const StrSet.Index) ListFmt {
    return .{
        .indices = indices,
        .str_set = str_set,
    };
}
