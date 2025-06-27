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
