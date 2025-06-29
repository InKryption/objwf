const std = @import("std");

pub const Indexer = std.ArrayHashMapUnmanaged(Index, void, HashCtx, true);

pub const Index = enum(Int) {
    pub const Int = u32;
    null = std.math.maxInt(Int),
    _,

    pub fn from(int: Int) Index {
        return @enumFromInt(int);
    }

    pub fn value(self: Index) ?Int {
        return switch (self) {
            .null => null,
            else => @intFromEnum(self),
        };
    }

    pub fn valueAllowNull(self: Index) Int {
        return @intFromEnum(self);
    }

    pub fn nonNull(self: Index) ?Index {
        return switch (self) {
            .null => null,
            else => |val| val,
        };
    }
};

const Hasher = std.hash.Wyhash;
const hasher_seed = 0;
fn hasherResultTruncate(result: u64) u32 {
    return @truncate(result);
}

pub fn hashCtx(str_bytes: []const u8) HashCtx {
    return .{ .str_bytes = str_bytes };
}
pub const HashCtx = struct {
    str_bytes: []const u8,

    pub fn hash(ctx: HashCtx, idx: Index) u32 {
        const str = idxSliceTo(ctx.str_bytes, idx);
        const adapter = hashStrAdapter(ctx.str_bytes);
        return adapter.hash(str);
    }

    pub fn eql(ctx: HashCtx, a: Index, b: Index, b_index: usize) bool {
        _ = ctx;
        _ = b_index;
        return a == b;
    }
};

pub fn hashStrAdapter(str_bytes: []const u8) HashStrAdapter {
    return .{ .str_bytes = str_bytes };
}
pub const HashStrAdapter = struct {
    str_bytes: []const u8,

    pub fn hash(ctx: HashStrAdapter, str: []const u8) u32 {
        _ = ctx;
        std.debug.assert(std.mem.indexOfScalar(u8, str, 0) == null);
        return hasherResultTruncate(Hasher.hash(hasher_seed, str));
    }

    pub fn eql(ctx: HashStrAdapter, a_str: []const u8, b: Index, b_index: usize) bool {
        _ = b_index;
        std.debug.assert(std.mem.indexOfScalar(u8, a_str, 0) == null);
        const b_str = idxSliceTo(ctx.str_bytes, b);
        return std.mem.eql(u8, a_str, b_str);
    }
};

pub fn hashIterAdapter(str_bytes: []const u8, comptime Iter: type) HashIterAdapter(Iter) {
    return .{ .str_bytes = str_bytes };
}
pub fn HashIterAdapter(comptime Iter: type) type {
    return struct {
        str_bytes: []const u8,
        const Self = @This();

        pub fn hash(_: Self, iter: Iter) u32 {
            var hasher: Hasher = .init(hasher_seed);

            iter.reset();
            while (iter.next()) |segment| {
                hasher.update(segment);
            }

            return hasherResultTruncate(hasher.final());
        }

        pub fn eql(ctx: Self, a_iter: Iter, b: Index, b_index: usize) bool {
            _ = b_index;
            const b_str = idxSliceTo(ctx.str_bytes, b);
            var index: usize = 0;

            a_iter.reset();
            while (a_iter.next()) |segment| {
                const str_rest = b_str[index..];
                if (segment.len > str_rest.len) return false;
                if (!std.mem.eql(u8, segment, str_rest[0..segment.len])) return false;
                index += segment.len;
            }
            std.debug.assert(index <= b_str.len);
            return index == b_str.len;
        }

        comptime {
            var err_msg: []const u8 = "";
            const IterNs = switch (@typeInfo(Iter)) {
                .pointer => |ptr_info| switch (ptr_info.size) {
                    .one => ptr_info.child,
                    else => Iter,
                },
                else => Iter,
            };

            if (!@hasDecl(IterNs, "reset")) {
                err_msg = err_msg ++ "- Missing `reset` method\n";
            } else switch (@TypeOf(Iter.reset)) {
                fn (Iter) void,
                fn (anytype) void,
                => {},
                else => |Reset| {
                    err_msg = err_msg ++ "- Expected `reset` declaration to be callable " ++
                        "as a method with no arguments returning void, " ++
                        "but it's actually `" ++ @typeName(Reset) ++ "`.";
                },
            }

            if (!@hasDecl(IterNs, "next")) {
                err_msg = err_msg ++ "- Missing `next` method\n";
            } else switch (@TypeOf(Iter.next)) {
                fn (Iter) ?[]const u8,
                fn (anytype) ?[]const u8,
                => {},
                else => |Next| {
                    err_msg = err_msg ++ "- Expected `reset` declaration to be callable " ++
                        "as a method with no arguments returning ?[]const u8, " ++
                        "but it's actually `" ++ @typeName(Next) ++ "`.";
                },
            }
            if (err_msg.len != 0) @compileError(
                "Problem(s) with Iter type `" ++
                    @typeName(Iter) ++ "`:\n" ++
                    err_msg,
            );
        }
    };
}

pub fn idxSliceTo(str_bytes: []const u8, idx: Index) []const u8 {
    const start = idx.value().?;
    return str_bytes[start..std.mem.indexOfScalarPos(u8, str_bytes, start, 0).?];
}
