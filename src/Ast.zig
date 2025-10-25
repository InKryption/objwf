const std = @import("std");
const Parser = @import("Parser.zig");
pub const StrSet = @import("StrSet.zig");

comptime {
    _ = Parser;
    _ = StrSet;
}

const Ast = @This();
nodes: std.MultiArrayList(Node).Slice,
extra: []const u32,
indices: Indices,
str_set: StrSet,

pub fn deinit(self: Ast, gpa: std.mem.Allocator) void {
    gpa.free(self.extra);
    self.str_set.deinit(gpa);
    self.indices.deinit(gpa);
    var nodes = self.nodes;
    nodes.deinit(gpa);
}

pub const empty: Ast = .{
    .indices = .empty,
    .extra = &.{},
    .str_set = .empty,
};

pub const ParseDiagnostic = Parser.ParseDiagnostic;
pub const ParseError = Parser.ParseError;

pub fn parse(
    gpa: std.mem.Allocator,
    src: []const u8,
    maybe_diag: ?*ParseDiagnostic,
) (ParseError || std.mem.Allocator.Error)!Ast {
    var nodes: std.MultiArrayList(Ast.Node) = .empty;
    defer nodes.deinit(gpa);

    var extra: std.ArrayListUnmanaged(u32) = .empty;
    defer extra.deinit(gpa);

    var indices_wip: Parser.IndicesWip = .empty;
    defer indices_wip.deinit(gpa);

    var str_set_wip: Parser.StrSetWip = .empty;
    defer str_set_wip.deinit(gpa);

    var parser: Parser = .{
        .lexer = .init(src),
        .line = 0,
        .debug = .init,
    };
    try parser.parse(gpa, maybe_diag, .{
        .nodes = &nodes,
        .extra = &extra,
        .indices_wip = &indices_wip,
        .str_set_wip = &str_set_wip,
    });

    var nodes_final = nodes.toOwnedSlice();
    errdefer nodes_final.deinit(gpa);

    const extra_final = try extra.toOwnedSlice(gpa);
    errdefer gpa.free(extra_final);

    const indices_final = try indices_wip.toOwnedIndices(gpa);
    errdefer indices_final.deinit(gpa);

    const str_set_final = try str_set_wip.toOwnedStrSet(gpa);
    errdefer str_set_final.deinit(gpa);

    return .{
        .nodes = nodes_final,
        .extra = extra_final,
        .indices = indices_final,
        .str_set = str_set_final,
    };
}

pub const Node = struct {
    tag: Tag,
    data: Data,

    pub const Index = u32;

    pub const Data = packed union {
        empty: enum (u32) { unset },
        index: u32,
        string: StrSet.Index,
        /// Source string for a value (ie an f32 string).
        value: StrSet.Index,
        ref_index: RefIndex,
        cstype: packed struct(u32) {
            value: CsType,
            _padding: enum(u28) { unset = 0, _ } = .unset,
        },
    };

    pub const Tag = enum(u8) {
        /// `call filename arg0 arg1` [...] `argN`
        /// filename = `str_set.getStr(.from(extra[data.index]))`
        /// arg[i] = `str_set.getStr(.from(extra[data.index + 1 + i]))` while `StrSet.Index.from(extra[data.index + 1 + i]) != .null`
        call,
        /// `call filename`
        /// filename = `str_set.getStr(data.string)`
        call_no_args,

        /// `mtllib filename`
        /// filename = `str_set.getStr(data.string)`
        mtllib_one,
        /// `mtllib filename0 filename1` [...] `filenameN`
        /// filename[i] = `str_set.getStr(.from(extra[data.index + 1 + i]))` while `StrSet.Index.from(extra[data.index + 1 + i]) != .null`
        mtllib_multi,

        /// `usemtl material_name`
        /// material_name = `str_set.getStr(data.string)`
        usemtl,

        /// `o name`
        /// name = `str_set.getStr(data.string)`
        o,

        /// `v x y z`
        /// x, y, z = `@as(*[3]StrSet.Index, @memCast(extra[data.index..][0..3])).*`
        /// The omitted `w` component is assumed to be 1.0
        v_xyz,
        /// `v x y z w`
        /// x, y, z, w = `@as(*[4]StrSet.Index, @memCast(extra[data.index..][0..4])).*`
        v_xyzw,

        /// `vn i j k`
        /// i, j, k = `@as(*[3]StrSet.Index, @memCast(extra[data.index..][0..3])).*`
        vn_ijk,

        /// `vp u`
        /// u = `str_set.getStr(data.value)`
        vp_u,
        /// `vp u v`
        /// u, v = `@as(*[2]StrSet.Index, @memCast(extra[data.index..][0..2])).*`
        vp_uv,
        /// `vp u v w`
        /// u, v, w = `@as(*[3]StrSet.Index, @memCast(extra[data.index..][0..3])).*`
        vp_uvw,

        /// `vn u`
        /// u = `str_set.getStr(data.value)`
        vt_u,
        /// `vn u v`
        /// u, v = `@as(*[2]StrSet.Index, @memCast(extra[data.index..][0..2])).*`
        vt_uv,
        /// `vn u v w`
        /// u, v, w = `@as(*[3]StrSet.Index, @memCast(extra[data.index..][0..3])).*`
        vt_uvw,

        /// `cstype rat type`
        /// rat = `data.cstype.value`
        cstype,

        /// `deg degu`
        /// degu = `str_set.getStr(data.value)`
        deg_curve,
        /// `deg degu degv`
        /// degu, degv = `@as(*[2]StrSet.Index, @memCast(extra[data.index..][0..2])).*`
        deg_surface,

        /// `g`
        g_empty,
        /// `g group`
        /// group = `str_set.getStr(data.string)`
        g_single,
        /// `g group0 group1` [...] `groupN`
        /// group[i] = `str_set.getStr(.from(extra[data.index + i]))` while `StrSet.Index.from(extra[data.index + i]) != .null`
        g_multi,

        /// `s group_number`
        /// group_number = `str_set.getStr(data.string)`
        s,

        /// `f v`
        /// v = `data.ref_index`
        f_one_v,
        /// `f v/vt`
        /// v, vt = `@as(RefIndex, @memCast(extra[data.index..][0..2])).*`
        f_one_v_vt,
        /// `f v//vn`
        /// v, vn = `@as(RefIndex, @memCast(extra[data.index..][0..2])).*`
        f_one_v_vn,
        /// `f v/vt/vn`
        /// v, vt, vn = `@as(RefIndex, @memCast(extra[data.index..][0..3])).*`
        f_one_v_vt_vn,

        /// `f v0 v1`
        /// v0, v1 = `@as(RefIndex, @memCast(extra[data.index..][0..2])).*`
        f_two_v,
        /// `f v0/vt0 v1/vt1`
        /// (v0, vt0), (v1, vt1) = `@as(RefIndex, @memCast(extra[data.index..][0..4])).*`
        f_two_v_vt,
        /// `f v0//vn0 v1//vn1`
        /// (v0, vn0), (v1, vn1) = `@as(RefIndex, @memCast(extra[data.index..][0..4])).*`
        f_two_v_vn,
        /// `f v0/vt0/vn0 v1/vt1/vn1`
        /// (v0, vt0, vn0), (v1, vt1, vn1) = `@as(RefIndex, @memCast(extra[data.index..][0..6])).*`
        f_two_v_vt_vn,

        /// `f v0 v1 v2`
        /// v0, v1, v2 = `@as(RefIndex, @memCast(extra[data.index..][0..3])).*`
        f_three_v,
        /// `f v0/vt0 v1/vt1 v2/vt2`
        /// (v0, vt0), (v1, vt1), (v2, vt2) = `@as(RefIndex, @memCast(extra[data.index..][0..6])).*`
        f_three_v_vt,
        /// `f v0//vn0 v1//vn1 v2//vn2`
        /// (v0, vn0), (v1, vn1), (v2, vn2) = `@as(RefIndex, @memCast(extra[data.index..][0..6])).*`
        f_three_v_vn,
        /// `f v0/vt0/vn0 v1/vt1/vn1 v2/vt2/vn2`
        /// (v0, vt0, vn0), (v1, vt1, vn1), (v2, vt2, vn2) = `@as(RefIndex, @memCast(extra[data.index..][0..9])).*`
        f_three_v_vt_vn,

        /// `f v0 v1 v2 v3`
        /// v0, v1, v2, v3 = `@as(RefIndex, @memCast(extra[data.index..][0..4])).*`
        f_four_v,
        /// `f v0/vt0 v1/vt1 v2/vt2 v3/vt3`
        /// v0, vt0, v1, vt1, v2, vt2, v3, vt3 = `@as(RefIndex, @memCast(extra[data.index..][0..8])).*`
        f_four_v_vt,
        /// `f v0//vn0 v1//vn1 v2//vn2 v3//vn3`
        /// (v0, vn0), (v1, vn1), (v2, vn2), (v3, vn3) = `@as(RefIndex, @memCast(extra[data.index..][0..8])).*`
        f_four_v_vn,
        /// `f v0/vt0/vn0 v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3`
        /// (v0, vt0, vn0), (v1, vt1, vn1), (v2, vt2, vn2), (v3, vt3, vn3) = `@as(RefIndex, @memCast(extra[data.index..][0..12])).*`
        f_four_v_vt_vn,

        /// `f v0 v1` [...] `vN`
        /// > while `RefIndex.from(extra[data.index + i]) != .null`
        /// v[i] = `RefIndex.from(extra[data.index + i])`
        f_multi_v,
        /// `f v0/vt0 v1/vt0` [...] `vN/vtN`
        /// > while `RefIndex.from(extra[data.index + i * 2 + 2]) != .null`
        /// v[i]  = `RefIndex.from(extra[data.index + i * 2 + 0])`
        /// vt[i] = `RefIndex.from(extra[data.index + i * 2 + 1])`
        f_multi_v_vt,
        /// `f v0//vn0 v1//vn0` [...] `vN//vnN`
        /// > while `RefIndex.from(extra[data.index + i * 2 + 2]) != .null`
        /// v [i] = `RefIndex.from(extra[data.index + i * 2 + 0])`
        /// vn[i] = `RefIndex.from(extra[data.index + i * 2 + 1])`
        f_multi_v_vn,
        /// `f v0/vt0/vn0 v1/vt1/vn0` [...] `vN/vtN/vnN`
        /// > while `RefIndex.from(extra[data.index + i * 2 + 0]) != .null`
        /// v[i]  = `RefIndex.from(extra[data.index + i * 2 + 0])`
        /// vt[i] = `RefIndex.from(extra[data.index + i * 2 + 1])`
        /// vn[i] = `RefIndex.from(extra[data.index + i * 2 + 2])`
        f_multi_v_vt_vn,
    };

    pub const CsType = packed struct(u4) {
        rational: bool,
        type: Type,

        pub const Type = enum(u3) {
            /// basis matrix
            bmatrix,
            /// Bezier
            bezier,
            /// B-spline
            bspline,
            /// Cardinal
            cardinal,
            /// Taylor
            taylor,
        };
    };

    pub const Call = struct {
        argv: []const StrSet.Index,

        /// May be null.
        pub fn filename(call: Call) StrSet.Index {
            return if (call.argv.len != 0) call.argv[0] else .null;
        }

        /// May be empty.
        pub fn args(call: Call) []const StrSet.Index {
            const start = @intFromBool(call.argv.len != 0);
            return call.argv[start..];
        }
    };

    pub const Unzipped = union(Tag) {
        call: Call,
        call_no_args: *const StrSet.Index,

        mtllib_one: *const StrSet.Index,
        mtllib_multi: []const StrSet.Index,

        usemtl: StrSet.Index,

        o: StrSet.Index,

        /// Like `v_xyzw`, but `w` defaults to `1.0`.
        v_xyz: [3]StrSet.Index,
        v_xyzw: [4]StrSet.Index,

        vn_ijk: [3]StrSet.Index,

        vp_u: StrSet.Index,
        vp_uv: [2]StrSet.Index,
        vp_uvw: [3]StrSet.Index,

        vt_u: StrSet.Index,
        vt_uv: [2]StrSet.Index,
        vt_uvw: [3]StrSet.Index,

        cstype: CsType,

        deg_curve: StrSet.Index,
        deg_surface: [2]StrSet.Index,

        g_empty,
        g_single: *const StrSet.Index,
        g_multi: []const StrSet.Index,

        s: StrSet.Index,

        f_one_v: *const RefIndex,
        f_one_v_vt: *const [2]RefIndex,
        f_one_v_vn: *const [2]RefIndex,
        f_one_v_vt_vn: *const [3]RefIndex,

        f_two_v: *const [2]RefIndex,
        f_two_v_vt: *const [2][2]RefIndex,
        f_two_v_vn: *const [2][2]RefIndex,
        f_two_v_vt_vn: *const [2][3]RefIndex,

        f_three_v: *const [3]RefIndex,
        f_three_v_vt: *const [3][2]RefIndex,
        f_three_v_vn: *const [3][2]RefIndex,
        f_three_v_vt_vn: *const [3][3]RefIndex,

        f_four_v: *const [4]RefIndex,
        f_four_v_vt: *const [4][2]RefIndex,
        f_four_v_vn: *const [4][2]RefIndex,
        f_four_v_vt_vn: *const [4][3]RefIndex,

        f_multi_v: []const RefIndex,
        f_multi_v_vt: []const [2]RefIndex,
        f_multi_v_vn: []const [2]RefIndex,
        f_multi_v_vt_vn: []const [3]RefIndex,

        pub fn full(self: Unzipped) Statement {
            return switch (self) {
                .call => |call| .{ .call = call },
                .call_no_args => |filename| .{ .call = .{ .argv = filename[0..1] } },

                .mtllib_one => |lib| .{ .mtllib = lib[0..1] },
                .mtllib_multi => |libs| .{ .mtllib = libs },
                .usemtl => |mtl| .{ .usemtl = mtl },

                .o => |name| .{ .o = name },

                .v_xyz => |xyz| .{ .v = xyz ++ .{.null} },
                .v_xyzw => |xyzw| .{ .v = xyzw },

                .vn_ijk => |ijk| .{ .vn = ijk },

                .vp_u => |u| .{ .vp = .{ u, .null, .null } },
                .vp_uv => |uv| .{ .vp = uv ++ .{.null} },
                .vp_uvw => |uvw| .{ .vp = uvw },

                .vt_u => |u| .{ .vt = .{ u, .null, .null } },
                .vt_uv => |uv| .{ .vt = uv ++ .{.null} },
                .vt_uvw => |uvw| .{ .vt = uvw },

                .cstype => |cstype| .{ .cstype = cstype },

                .deg_curve => |degu| .{ .deg = .{ degu, .null } },
                .deg_surface => |deguv| .{ .deg = deguv },

                .g_empty => .{ .g = &.{} },
                .g_single => |group| .{ .g = group[0..1] },
                .g_multi => |groups| .{ .g = groups },

                .s => |group_number| .{ .s = group_number },

                .f_one_v => |one_v| .{ .f = .{ .vertex_only = @as(*const RefIndex, one_v)[0..1] } },
                .f_one_v_vt => |one_v_vt| .{ .f = .{ .only_vt = @ptrCast(one_v_vt) } },
                .f_one_v_vn => |one_v_vn| .{ .f = .{ .only_vn = @ptrCast(one_v_vn) } },
                .f_one_v_vt_vn => |one_v_vt_vn| .{ .f = .{ .both_vt_vn = @ptrCast(one_v_vt_vn) } },

                .f_two_v => |two_v| .{ .f = .{ .vertex_only = two_v } },
                .f_two_v_vt => |two_v_vt| .{ .f = .{ .only_vt = two_v_vt } },
                .f_two_v_vn => |two_v_vn| .{ .f = .{ .only_vn = two_v_vn } },
                .f_two_v_vt_vn => |two_v_vt_vn| .{ .f = .{ .both_vt_vn = two_v_vt_vn } },

                .f_three_v => |three_v| .{ .f = .{ .vertex_only = three_v } },
                .f_three_v_vt => |three_v_vt| .{ .f = .{ .only_vt = three_v_vt } },
                .f_three_v_vn => |three_v_vn| .{ .f = .{ .only_vn = three_v_vn } },
                .f_three_v_vt_vn => |three_v_vt_vn| .{ .f = .{ .both_vt_vn = three_v_vt_vn } },

                .f_four_v => |four_v| .{ .f = .{ .vertex_only = four_v } },
                .f_four_v_vt => |four_v_vt| .{ .f = .{ .only_vt = four_v_vt } },
                .f_four_v_vn => |four_v_vn| .{ .f = .{ .only_vn = four_v_vn } },
                .f_four_v_vt_vn => |four_v_vt_vn| .{ .f = .{ .both_vt_vn = four_v_vt_vn } },

                .f_multi_v => |multi_v| .{ .f = .{ .vertex_only = multi_v } },
                .f_multi_v_vt => |multi_v_vt| .{ .f = .{ .only_vt = multi_v_vt } },
                .f_multi_v_vn => |multi_v_vn| .{ .f = .{ .only_vn = multi_v_vn } },
                .f_multi_v_vt_vn => |multi_v_vt_vn| .{ .f = .{ .both_vt_vn = multi_v_vt_vn } },
            };
        }
    };
};

pub const Statement = union(Kind) {
    call: Node.Call,
    mtllib: []const StrSet.Index,
    usemtl: StrSet.Index,
    /// Never null.
    o: StrSet.Index,
    /// x, y, z, w = v
    /// x, y and z are never null.
    /// w may be null, in which case it should be interpreted as `1.0`.
    v: [4]StrSet.Index,
    /// i, j, k = vp
    /// i, j, and k are never null.
    vn: [3]StrSet.Index,
    /// u, v, w = vp
    /// u is never null.
    /// v may be null iff w is null.
    /// w may be null.
    vp: [3]StrSet.Index,
    /// u, v, w = vp
    /// u is never null.
    /// v may be null iff w is null.
    /// w may be null.
    vt: [3]StrSet.Index,
    cstype: Node.CsType,
    /// degu, degv = deg
    /// degu is never null.
    /// degv may be null.
    ///
    /// if degv is null, it's a curve.
    /// if degv is not null, it's a surface.
    deg: [2]StrSet.Index,
    g: []const StrSet.Index,
    s: StrSet.Index,
    f: FaceTriplet.List,

    pub const Kind = enum {
        call,
        mtllib,
        usemtl,
        o,
        v,
        vn,
        vp,
        vt,
        cstype,
        deg,
        g,
        s,
        f,
    };
};

/// Reference to an indexed vertex data.
/// 0-based reference to enumerated vertex data.
pub const RefIndex = enum(Int) {
    pub const Int = u32;
    null = std.math.maxInt(Int),
    _,

    /// Asserts `self != .null`, and returns `self`.
    pub fn nonNull(self: RefIndex) ?RefIndex {
        return switch (self) {
            .null => null,
            else => |val| val,
        };
    }

    pub fn from(int: Int) RefIndex {
        return @enumFromInt(int);
    }

    pub const OneBased = std.math.IntFittingRange(
        1,
        std.math.maxInt(Int) + 1,
    );
    /// Asserts `int != 0`.
    /// Asserts `int <= std.math.maxInt(Int) + 1`.
    pub fn fromOneBased(int: OneBased) RefIndex {
        std.debug.assert(int != 0);
        const casted: Int = @intCast(int - 1);
        return @enumFromInt(casted);
    }

    pub fn value(self: RefIndex) ?Int {
        return switch (self) {
            .null => null,
            else => @intFromEnum(self),
        };
    }

    pub fn valueAllowNull(self: RefIndex) Int {
        return @intFromEnum(self);
    }
};

pub const FaceTriplet = struct {
    /// Cannot be `.null`.
    v: RefIndex,
    vt: RefIndex,
    vn: RefIndex,

    pub const Layout = enum {
        /// `v`
        vertex_only,
        /// `v/vt`
        only_vt,
        /// `v//vn`
        only_vn,
        /// `v/vt/vn`
        both_vt_vn,

        pub fn zippedLength(self: Layout) u2 {
            return switch (self) {
                .vertex_only => 1,
                .only_vt => 2,
                .only_vn => 2,
                .both_vt_vn => 3,
            };
        }
    };

    pub fn layout(self: FaceTriplet) Layout {
        std.debug.assert(self.v != .null);
        if (self.vt == .null and self.vn == .null) return .vertex_only;
        if (self.vt == .null) return .only_vn;
        if (self.vn == .null) return .only_vt;
        return .both_vt_vn;
    }

    pub const List = union(Layout) {
        /// `v`
        vertex_only: []const RefIndex,
        /// `v/vt`
        only_vt: []const [2]RefIndex,
        /// `v//vn`
        only_vn: []const [2]RefIndex,
        /// `v/vt/vn`
        both_vt_vn: []const [3]RefIndex,

        pub fn fmt(list: List) Fmt {
            return .{ .list = list };
        }

        pub const Fmt = struct {
            list: List,

            pub fn format(
                self: Fmt,
                writer: *std.Io.Writer,
            ) std.io.Writer.Error!void {
                switch (self.list) {
                    inline else => |indices, tag| for (indices, 0..) |index_or_list, triplet_i| {
                        if (triplet_i != 0) try writer.writeAll(" ");
                        switch (tag) {
                            .vertex_only => try writer.printInt(index_or_list.value().?, 10, .lower, .{}),
                            .only_vt => for (index_or_list, 0..) |ref_index, ref_index_i| {
                                if (ref_index_i != 0) try writer.writeByte("/");
                                try writer.printInt(ref_index.value().?, 10, .lower, .{});
                            },
                            .only_vn => for (index_or_list, 0..) |ref_index, ref_index_i| {
                                if (ref_index_i != 0) try writer.writeAll("//");
                                try writer.printInt(ref_index.value().?, 10, .lower, .{});
                            },
                            .both_vt_vn => for (index_or_list, 0..) |ref_index, ref_index_i| {
                                if (ref_index_i != 0) try writer.writeAll("/");
                                try writer.printInt(ref_index.value().?, 10, .lower, .{});
                            },
                        }
                    },
                }
            }
        };
    };
};

pub const Indices = struct {
    list: []const Node.Index,
    offsets: Offsets,

    pub const empty: Indices = .{
        .list = &.{},
        .offsets = .empty,
    };

    pub fn deinit(self: Indices, gpa: std.mem.Allocator) void {
        gpa.free(self.list);
    }

    pub fn endOffsetOf(self: Indices, comptime field: Offsets.FieldTag) u32 {
        if (comptime Offsets.endFieldOf(field)) |end_tag| {
            return @field(self.offsets, @tagName(end_tag));
        } else {
            return @intCast(self.list.len);
        }
    }

    pub const Offsets = struct {
        /// Vertex Positions
        comptime positions: u32 = 0,
        /// Texture Coordinates
        textures: u32,
        /// Vertex Normals
        normals: u32,
        /// Parameter Space Vertices
        parameter_spaces: u32,

        pub const empty: Offsets = .{
            .positions = 0,
            .textures = 0,
            .normals = 0,
            .parameter_spaces = 0,
        };

        pub const FieldTag = std.meta.FieldEnum(Offsets);

        /// Returns the field tag such that `list[offsets.<start_field>..offsets.<end_field>]` is
        /// the list of the indices corresponding to the name of `start_field`.
        /// Returns null for the last field, where the list is `list[offsets.<start_field>..]`.
        pub fn endFieldOf(start_field: FieldTag) ?FieldTag {
            return switch (start_field) {
                .positions => //
                .textures,
                .textures => //
                .normals,
                .normals => //
                .parameter_spaces,
                .parameter_spaces => //
                null,
            };
        }
    };

    pub fn slices(self: Indices) Slices {
        const offsets = self.offsets;
        return .{
            // zig fmt: off
            .positions        = self.list[offsets.positions       ..offsets.textures],
            .textures         = self.list[offsets.textures        ..offsets.normals],
            .normals          = self.list[offsets.normals         ..offsets.parameter_spaces],
            .parameter_spaces = self.list[offsets.parameter_spaces..],
            // zig fmt: on
        };
    }

    pub const Slices = struct {
        /// Vertex Positions
        positions: []const Node.Index,
        /// Texture Coordinates
        textures: []const Node.Index,
        /// Vertex Normals
        normals: []const Node.Index,
        /// Parameter Space Vertices
        parameter_spaces: []const Node.Index,

        pub const empty: Slices = .{
            .positions = &.{},
            .textures = &.{},
            .normals = &.{},
            .parameter_spaces = &.{},
        };
    };
};

pub fn unzipped(ast: Ast, node_index: Node.Index) Node.Unzipped {
    const node_tag: *const Node.Tag = &ast.nodes.items(.tag)[node_index];
    const node_data: *const Node.Data = &ast.nodes.items(.data)[node_index];
    return switch (node_tag.*) {
        .g_empty,
        => .g_empty,

        .cstype,
        => .{ .cstype = node_data.cstype.value },

        .g_multi,
        => .{ .g_multi = ast.strListSliceToNull(node_data.index) },

        .call,
        => .{ .call = .{
            .argv = ast.strListSliceToNull(node_data.index),
        } },

        .mtllib_one => .{ .mtllib_one = &node_data.string },
        .mtllib_multi => .{ .mtllib_multi = ast.strListSliceToNull(node_data.index) },
        .usemtl => .{ .usemtl = node_data.string },

        inline //
        .vp_u,
        .vt_u,
        .deg_curve,
        => |tag| @unionInit(Node.Unzipped, @tagName(tag), node_data.value),

        inline //
        .o,
        .s,
        => |tag| @unionInit(Node.Unzipped, @tagName(tag), node_data.string),
        inline //
        .call_no_args,
        .g_single,
        => |tag| @unionInit(Node.Unzipped, @tagName(tag), &node_data.string),

        inline //
        .deg_surface,

        .v_xyz,
        .v_xyzw,

        .vp_uv,
        .vp_uvw,

        .vt_uv,
        .vt_uvw,

        .vn_ijk,
        => |tag| blk: {
            const len = comptime switch (tag) {
                .deg_surface => 2,

                .v_xyz => 3,
                .v_xyzw => 4,

                .vp_uv => 2,
                .vp_uvw => 3,

                .vt_uv => 2,
                .vt_uvw => 3,

                .vn_ijk => 3,
                else => unreachable,
            };
            const data: *const [len]StrSet.Index =
                elemCast(StrSet.Index, ast.extra[node_data.index..][0..len]);
            break :blk @unionInit(Node.Unzipped, @tagName(tag), data.*);
        },

        inline //
        .f_one_v,
        .f_one_v_vt,
        .f_one_v_vn,
        .f_one_v_vt_vn,

        .f_two_v,
        .f_two_v_vt,
        .f_two_v_vn,
        .f_two_v_vt_vn,

        .f_three_v,
        .f_three_v_vt,
        .f_three_v_vn,
        .f_three_v_vt_vn,

        .f_four_v,
        .f_four_v_vt,
        .f_four_v_vn,
        .f_four_v_vt_vn,
        => |tag| blk: {
            const len, const elem_size = comptime switch (tag) {
                // zig fmt: off
                    .f_one_v         => .{ 1, 1 },
                    .f_one_v_vt      => .{ 1, 2 },
                    .f_one_v_vn      => .{ 1, 2 },
                    .f_one_v_vt_vn   => .{ 1, 3 },

                    .f_two_v         => .{ 2, 1 },
                    .f_two_v_vt      => .{ 2, 2 },
                    .f_two_v_vn      => .{ 2, 2 },
                    .f_two_v_vt_vn   => .{ 2, 3 },

                    .f_three_v       => .{ 3, 1 },
                    .f_three_v_vt    => .{ 3, 2 },
                    .f_three_v_vn    => .{ 3, 2 },
                    .f_three_v_vt_vn => .{ 3, 3 },

                    .f_four_v        => .{ 4, 1 },
                    .f_four_v_vt     => .{ 4, 2 },
                    .f_four_v_vn     => .{ 4, 2 },
                    .f_four_v_vt_vn  => .{ 4, 3 },
                    else => unreachable,
                    // zig fmt: on
            };
            const Elem = if (elem_size == 1) RefIndex else [elem_size]RefIndex;
            const data: *const [len]Elem = elemCast(Elem, ast.extra[node_data.index..][0..len]);
            break :blk @unionInit(Node.Unzipped, @tagName(tag), if (len == 1) &data[0] else data);
        },

        inline //
        .f_multi_v,
        .f_multi_v_vt,
        .f_multi_v_vn,
        .f_multi_v_vt_vn,
        => |tag| blk: {
            const layout: FaceTriplet.Layout = comptime switch (tag) {
                .f_multi_v => .vertex_only,
                .f_multi_v_vt => .only_vt,
                .f_multi_v_vn => .only_vn,
                .f_multi_v_vt_vn => .both_vt_vn,
                else => unreachable,
            };
            break :blk @unionInit(
                Node.Unzipped,
                @tagName(tag),
                ast.multiFaceList(node_data.index, layout),
            );
        },
    };
}

pub fn full(self: Ast, node_index: Node.Index) Statement {
    return self.unzipped(node_index).full();
}

fn strListSliceToNull(ast: Ast, index: u32) []const Ast.StrSet.Index {
    return std.mem.sliceTo(elemCast(Ast.StrSet.Index, ast.extra[index..]), .null);
}

fn refIdxListSliceToNull(ast: Ast, index: u32) []const Ast.RefIndex {
    return std.mem.sliceTo(elemCast(Ast.RefIndex, ast.extra[index..]), .null);
}

fn multiFaceList(
    self: Ast,
    index: u32,
    comptime layout: FaceTriplet.Layout,
) []const if (layout.zippedLength() != 1) [layout.zippedLength()]RefIndex else RefIndex {
    const raw: []const RefIndex = self.refIdxListSliceToNull(index);
    const elem_len = comptime layout.zippedLength();
    std.debug.assert(raw.len % elem_len == 0);
    const Elem = if (elem_len == 1) RefIndex else [elem_len]RefIndex;
    return elemCast(Elem, raw);
}

/// TODO: replace with @memCast or w/e when that gets added
fn elemCast(comptime Elem: type, src: anytype) @Type(.{ .pointer = blk: {
    const Src = @TypeOf(src);
    var info = @typeInfo(Src).pointer;
    switch (info.size) {
        .one => switch (@typeInfo(info.child)) {
            .array => |array_info| {
                info.child = [array_info.len]Elem;
            },
            else => info.child = Elem,
        },
        .many, .slice => info.child = Elem,
        .c => @compileError("C pointers unsupported"),
    }
    break :blk info;
} }) {
    return @ptrCast(src);
}

const TestUnzippedNode = union(Statement.Kind) {
    call: struct { ?[]const u8, []const []const u8 },
    mtllib: []const []const u8,
    usemtl: ?[]const u8,

    o: ?[]const u8,

    v: struct { []const u8, []const u8, []const u8, ?[]const u8 },
    vn: [3][]const u8,
    vp: struct { []const u8, ?[]const u8, ?[]const u8 },
    vt: struct { []const u8, ?[]const u8, ?[]const u8 },
    cstype: Node.CsType,
    deg: struct { []const u8, ?[]const u8 },
    g: []const []const u8,
    s: ?[]const u8,
    f: FaceTriplet.List,
};

fn testAstParse(
    src: []const u8,
    expected: []const TestUnzippedNode,
) !void {
    const gpa = std.testing.allocator;

    var diag: ParseDiagnostic = undefined;
    const ast = Ast.parse(gpa, src, &diag) catch |err| return switch (err) {
        error.OutOfMemory => |e| e,
        else => |e| {
            std.log.err("{f}", .{diag.tagged(e).fmt(src)});
        },
    };
    defer ast.deinit(gpa);

    var actual: std.ArrayListUnmanaged(TestUnzippedNode) = .empty;
    defer actual.deinit(gpa);
    try actual.ensureTotalCapacityPrecise(gpa, expected.len);

    var node_arena_state: std.heap.ArenaAllocator = .init(gpa);
    defer node_arena_state.deinit();
    const node_arena = node_arena_state.allocator();

    for (0..ast.nodes.len) |node_index| {
        try actual.append(gpa, switch (ast.full(@intCast(node_index))) {
            inline else => |payload, tag| @unionInit(TestUnzippedNode, @tagName(tag), switch (tag) {
                .call => .{
                    ast.str_set.getStr(payload.filename()),
                    try strIndicesToStrings(node_arena, ast.str_set, payload.args()),
                },
                .mtllib => try strIndicesToStrings(node_arena, ast.str_set, payload),
                .usemtl => ast.str_set.getStr(payload),
                .o => ast.str_set.getStr(payload),
                .v => v: {
                    const x: StrSet.Index, //
                    const y: StrSet.Index, //
                    const z: StrSet.Index, //
                    const w: StrSet.Index //
                    = payload;
                    break :v .{
                        ast.str_set.getStr(x).?,
                        ast.str_set.getStr(y).?,
                        ast.str_set.getStr(z).?,
                        ast.str_set.getStr(w),
                    };
                },
                .vn => vn: {
                    const i: StrSet.Index, //
                    const j: StrSet.Index, //
                    const k: StrSet.Index //
                    = payload;
                    break :vn .{
                        ast.str_set.getStr(i).?,
                        ast.str_set.getStr(j).?,
                        ast.str_set.getStr(k).?,
                    };
                },
                .vp => vp: {
                    const u: StrSet.Index, //
                    const v: StrSet.Index, //
                    const w: StrSet.Index //
                    = payload;
                    break :vp .{
                        ast.str_set.getStr(u).?,
                        ast.str_set.getStr(v),
                        ast.str_set.getStr(w),
                    };
                },
                .vt => vt: {
                    const u: StrSet.Index, //
                    const v: StrSet.Index, //
                    const w: StrSet.Index //
                    = payload;
                    break :vt .{
                        ast.str_set.getStr(u).?,
                        ast.str_set.getStr(v),
                        ast.str_set.getStr(w),
                    };
                },
                .cstype => payload,
                .deg => deg: {
                    const degu: StrSet.Index, const degv: StrSet.Index = payload;
                    break :deg .{
                        ast.str_set.getStr(degu).?,
                        ast.str_set.getStr(degv),
                    };
                },
                .g => try strIndicesToStrings(node_arena, ast.str_set, payload),

                .s => ast.str_set.getStr(payload),
                .f => payload,
            }),
        });
    }

    try std.testing.expectEqualDeep(expected, actual.items);
}

const TestDiagnostic = union(ParseDiagnostic.Tagged.Tag) {
    InvalidBytes: LineSrc,
    InvalidEscape: LineSrc,
    InvalidParameter: LineSrc,
    InvalidPrefix: LineSrc,
    InvalidReferenceIndexTriplet: LineSrc,
    IncompleteStatement: LineSrc,
    UnexpectedParameter: LineSrc,
    MixedFaceLayouts: MixedFaceLayoutsPayload,

    const LineSrc = struct {
        line: u32,
        src: []const u8,
    };

    const MixedFaceLayoutsPayload = struct {
        first: Face,
        bad: Face,

        const Face = struct {
            /// line
            u32,
            /// layout
            FaceTriplet.Layout,
            /// src
            []const u8,
        };
    };
};

fn testAstDiagnostic(
    src: []const u8,
    expected: TestDiagnostic,
) !void {
    var diag: ParseDiagnostic = undefined;
    if (Ast.parse(std.testing.allocator, src, &diag)) |ast| {
        defer ast.deinit(std.testing.allocator);
        try std.testing.expectEqual(expected, null);
    } else |parse_err| switch (parse_err) {
        error.OutOfMemory => |oom| return oom,
        else => |diag_err| {
            if (@errorReturnTrace()) |ert| {
                ert.index = 0;
            }
            const actual: TestDiagnostic = switch (diag.tagged(diag_err)) {
                inline //
                .InvalidBytes,
                .InvalidEscape,
                .InvalidParameter,
                .InvalidPrefix,
                .InvalidReferenceIndexTriplet,
                .IncompleteStatement,
                .UnexpectedParameter,
                => |line_loc, tag| @unionInit(TestDiagnostic, @tagName(tag), .{
                    .line = line_loc.line,
                    .src = line_loc.loc.getStr(src),
                }),

                .MixedFaceLayouts => |mfl| .{ .MixedFaceLayouts = .{
                    .first = .{ mfl.first.line, mfl.first.layout, mfl.first.loc.getStr(src) },
                    .bad = .{ mfl.bad.line, mfl.bad.layout, mfl.bad.loc.getStr(src) },
                } },
            };
            try std.testing.expectEqual(
                std.meta.activeTag(expected),
                std.meta.activeTag(actual),
            );
            switch (expected) {
                inline else => |expected_pl, tag| {
                    const actual_pl = @field(actual, @tagName(tag));
                    switch (@TypeOf(expected_pl)) {
                        TestDiagnostic.LineSrc => {
                            try std.testing.expectEqual(expected_pl.line, actual_pl.line);
                            try std.testing.expectEqualStrings(expected_pl.src, actual_pl.src);
                        },
                        TestDiagnostic.MixedFaceLayoutsPayload => {
                            const expected_first_line: u32, //
                            const expected_first_layout: FaceTriplet.Layout, //
                            const expected_first_src: []const u8 //
                            = expected_pl.first;
                            const expected_bad_line: u32, //
                            const expected_bad_layout: FaceTriplet.Layout, //
                            const expected_bad_src: []const u8 //
                            = expected_pl.bad;

                            const actual_first_line: u32, //
                            const actual_first_layout: FaceTriplet.Layout, //
                            const actual_first_src: []const u8 //
                            = actual_pl.first;
                            const actual_bad_line: u32, //
                            const actual_bad_layout: FaceTriplet.Layout, //
                            const actual_bad_src: []const u8 //
                            = actual_pl.bad;

                            var err: error{TestExpectedEqual}!void = {};
                            std.testing.expectEqual(expected_first_line, actual_first_line) catch |e| {
                                err = e;
                            };
                            std.testing.expectEqual(expected_first_layout, actual_first_layout) catch |e| {
                                err = e;
                            };
                            std.testing.expectEqualStrings(expected_first_src, actual_first_src) catch |e| {
                                err = e;
                            };

                            std.testing.expectEqual(expected_bad_line, actual_bad_line) catch |e| {
                                err = e;
                            };
                            std.testing.expectEqual(expected_bad_layout, actual_bad_layout) catch |e| {
                                err = e;
                            };
                            std.testing.expectEqualStrings(expected_bad_src, actual_bad_src) catch |e| {
                                err = e;
                            };

                            return err;
                        },
                        else => comptime unreachable,
                    }
                },
            }
        },
    }
}

fn strIndicesToStrings(
    gpa: std.mem.Allocator,
    str_set: StrSet,
    indices: []const StrSet.Index,
) std.mem.Allocator.Error![]const []const u8 {
    var list: std.ArrayListUnmanaged([]const u8) = .empty;
    defer list.deinit(gpa);
    try list.ensureTotalCapacityPrecise(gpa, indices.len);
    std.debug.assert(list.unusedCapacitySlice().len == 0);
    for (indices) |str_idx| list.appendAssumeCapacity(str_set.getStr(str_idx).?);
    return try list.toOwnedSlice(gpa);
}

test "Empty Obj whitespace" {
    try testAstParse("", &.{});
    try testAstParse("\n", &.{});
    try testAstParse("\n" ** 23, &.{});
}

test "Empty Obj with comments" {
    try testAstParse(
        \\# foo bar baz
        \\# look      \
        \\  at        \
        \\  this      \
        \\  multiline \
        \\  comment   \
        \\  using     \
        \\  escaped   \
        \\  newlines
        \\
        \\#fizz buzz
        \\
    ,
        &.{},
    );
}

test "Escaped newlines" {
    try testAstParse( // chained newline escapes
        \\d\
        \\e\
        \\g\
        \\ 2\
        \\\
        \\.0
        \\o\
        \\ \
        \\foo\
        \\bar\
        \\\
        \\\
        \\\
        \\baz
    ,
        &.{
            .{ .deg = .{ "2.0", null } },
            .{ .o = "foobarbaz" },
        },
    );
}

test "Vertices" {
    try testAstParse(
        \\v \
        \\  -5.000000 \
        \\   5.000000 \
        \\   0.000000
    ,
        &.{
            .{ .v = .{ "-5.000000", "5.000000", "0.000000", null } },
        },
    );

    try testAstParse(
        \\# Polygonal and free-form geometry statements.
        \\v \
        \\  -5.000000 \
        \\   5.000000 \
        \\   0.000000
        \\v \
        \\  -5.000000 \
        \\  -5.000000 \
        \\   0.000000
        \\v \
        \\   5.000000 \
        \\  -5.000000 \
        \\   0.000000
        \\v \
        \\   5.000000 \
        \\   5.000000 \
        \\   0.000000
        \\
        \\# Vertex statements for both polygonal and free-form geometry.
        \\vt     -5.000000       5.000000       0.000000
        \\vt     -5.000000      -5.000000       0.000000
        \\vt      5.000000      -5.000000       0.000000
        \\vt      5.000000       5.000000       0.000000
        \\
        \\# Vertex statements for polygonal and free-form geometry.
        \\vn      0.000000       0.000000       1.000000
        \\vn      0.000000       0.000000       1.000000
        \\vn      0.000000       0.000000       1.000000
        \\vn      0.000000       0.000000       1.000000
        \\
        \\# Vertex statements for free-form geometry.
        \\vp      0.210000       3.590000
        \\vp      0.000000       0.000000
        \\vp      1.000000       0.000000
        \\vp      0.500000       0.500000
    ,
        &.{
            .{ .v = .{ "-5.000000", "5.000000", "0.000000", null } },
            .{ .v = .{ "-5.000000", "-5.000000", "0.000000", null } },
            .{ .v = .{ "5.000000", "-5.000000", "0.000000", null } },
            .{ .v = .{ "5.000000", "5.000000", "0.000000", null } },

            .{ .vt = .{ "-5.000000", "5.000000", "0.000000" } },
            .{ .vt = .{ "-5.000000", "-5.000000", "0.000000" } },
            .{ .vt = .{ "5.000000", "-5.000000", "0.000000" } },
            .{ .vt = .{ "5.000000", "5.000000", "0.000000" } },

            .{ .vn = .{ "0.000000", "0.000000", "1.000000" } },
            .{ .vn = .{ "0.000000", "0.000000", "1.000000" } },
            .{ .vn = .{ "0.000000", "0.000000", "1.000000" } },
            .{ .vn = .{ "0.000000", "0.000000", "1.000000" } },

            .{ .vp = .{ "0.210000", "3.590000", null } },
            .{ .vp = .{ "0.000000", "0.000000", null } },
            .{ .vp = .{ "1.000000", "0.000000", null } },
            .{ .vp = .{ "0.500000", "0.500000", null } },
        },
    );
}

test "curves/surfaces" {
    try testAstParse(
        \\deg 5.0
        \\deg 5.0 5.0
        \\deg 5.0 3.0
        \\deg 3.0 5.0
        \\deg 3.0 3.0
        \\deg 3.0
    ,
        &.{
            .{ .deg = .{ "5.0", null } },
            .{ .deg = .{ "5.0", "5.0" } },
            .{ .deg = .{ "5.0", "3.0" } },
            .{ .deg = .{ "3.0", "5.0" } },
            .{ .deg = .{ "3.0", "3.0" } },
            .{ .deg = .{ "3.0", null } },
        },
    );
}

test "cstype" {
    try testAstParse(
        \\cstype rat bmatrix
        \\cstype taylor
    ,
        &.{
            .{ .cstype = .{ .rational = true, .type = .bmatrix } },
            .{ .cstype = .{ .rational = false, .type = .taylor } },
        },
    );
}

test "faces" {
    // -- one -- //
    try testAstParse(
        \\f 1
        \\f 1/2
        \\f 1//2
        \\f 1/2/3
        \\
    ,
        &.{
            .{ .f = .{ .vertex_only = &.{.fromOneBased(1)} } },
            .{ .f = .{ .only_vt = &.{.{ .fromOneBased(1), .fromOneBased(2) }} } },
            .{ .f = .{ .only_vn = &.{.{ .fromOneBased(1), .fromOneBased(2) }} } },
            .{ .f = .{ .both_vt_vn = &.{.{ .fromOneBased(1), .fromOneBased(2), .fromOneBased(3) }} } },
        },
    );

    // -- two -- //
    try testAstParse(
        \\f 1 2
        \\f 1/2 3/4
        \\f 1//2 3//4
        \\f 1/2/3 4/5/6
        \\
    ,
        &.{
            .{ .f = .{ .vertex_only = &.{
                .fromOneBased(1),
                .fromOneBased(2),
            } } },
            .{ .f = .{ .only_vt = &.{
                .{ .fromOneBased(1), .fromOneBased(2) },
                .{ .fromOneBased(3), .fromOneBased(4) },
            } } },
            .{ .f = .{ .only_vn = &.{
                .{ .fromOneBased(1), .fromOneBased(2) },
                .{ .fromOneBased(3), .fromOneBased(4) },
            } } },
            .{ .f = .{ .both_vt_vn = &.{
                .{ .fromOneBased(1), .fromOneBased(2), .fromOneBased(3) },
                .{ .fromOneBased(4), .fromOneBased(5), .fromOneBased(6) },
            } } },
        },
    );

    // -- three -- //
    try testAstParse(
        \\f 1 2 3
        \\f 1/2 3/4 5/6
        \\f 1//2 3//4 5//6
        \\f 1/2/3 4/5/6 7/8/9
        \\
    ,
        &.{
            .{ .f = .{ .vertex_only = &.{
                .fromOneBased(1),
                .fromOneBased(2),
                .fromOneBased(3),
            } } },
            .{ .f = .{ .only_vt = &.{
                .{ .fromOneBased(1), .fromOneBased(2) },
                .{ .fromOneBased(3), .fromOneBased(4) },
                .{ .fromOneBased(5), .fromOneBased(6) },
            } } },
            .{ .f = .{ .only_vn = &.{
                .{ .fromOneBased(1), .fromOneBased(2) },
                .{ .fromOneBased(3), .fromOneBased(4) },
                .{ .fromOneBased(5), .fromOneBased(6) },
            } } },
            .{ .f = .{ .both_vt_vn = &.{
                .{ .fromOneBased(1), .fromOneBased(2), .fromOneBased(3) },
                .{ .fromOneBased(4), .fromOneBased(5), .fromOneBased(6) },
                .{ .fromOneBased(7), .fromOneBased(8), .fromOneBased(9) },
            } } },
        },
    );

    // -- four -- //
    try testAstParse(
        "f 1 2 3 4",
        &.{.{ .f = .{ .vertex_only = &.{
            .fromOneBased(1),
            .fromOneBased(2),
            .fromOneBased(3),
            .fromOneBased(4),
        } } }},
    );
    try testAstParse(
        "f 1/2 3/4 5/6 7/8",
        &.{.{ .f = .{ .only_vt = &.{
            .{ .fromOneBased(1), .fromOneBased(2) },
            .{ .fromOneBased(3), .fromOneBased(4) },
            .{ .fromOneBased(5), .fromOneBased(6) },
            .{ .fromOneBased(7), .fromOneBased(8) },
        } } }},
    );
    try testAstParse(
        "f 1//2 3//4 5//6 7//8",
        &.{.{ .f = .{ .only_vn = &.{
            .{ .fromOneBased(1), .fromOneBased(2) },
            .{ .fromOneBased(3), .fromOneBased(4) },
            .{ .fromOneBased(5), .fromOneBased(6) },
            .{ .fromOneBased(7), .fromOneBased(8) },
        } } }},
    );
    try testAstParse(
        "f 1/2/3 4/5/6 7/8/9 10/11/12",
        &.{.{ .f = .{ .both_vt_vn = &.{
            .{ .fromOneBased(1), .fromOneBased(2), .fromOneBased(3) },
            .{ .fromOneBased(4), .fromOneBased(5), .fromOneBased(6) },
            .{ .fromOneBased(7), .fromOneBased(8), .fromOneBased(9) },
            .{ .fromOneBased(10), .fromOneBased(11), .fromOneBased(12) },
        } } }},
    );

    // -- five -- //
    try testAstParse("f 1 2 3 4 5", &.{.{ .f = .{ .vertex_only = &.{
        .fromOneBased(1),
        .fromOneBased(2),
        .fromOneBased(3),
        .fromOneBased(4),
        .fromOneBased(5),
    } } }});
    try testAstParse("f 1/2 3/4 5/6 7/8 9/10", &.{.{ .f = .{ .only_vt = &.{
        .{ .fromOneBased(1), .fromOneBased(2) },
        .{ .fromOneBased(3), .fromOneBased(4) },
        .{ .fromOneBased(5), .fromOneBased(6) },
        .{ .fromOneBased(7), .fromOneBased(8) },
        .{ .fromOneBased(9), .fromOneBased(10) },
    } } }});
    try testAstParse("f 1//2 3//4 5//6 7//8 9//10", &.{.{ .f = .{ .only_vn = &.{
        .{ .fromOneBased(1), .fromOneBased(2) },
        .{ .fromOneBased(3), .fromOneBased(4) },
        .{ .fromOneBased(5), .fromOneBased(6) },
        .{ .fromOneBased(7), .fromOneBased(8) },
        .{ .fromOneBased(9), .fromOneBased(10) },
    } } }});
    try testAstParse("f 1/2/3 4/5/6 7/8/9 10/11/12 13/14/15", &.{.{ .f = .{ .both_vt_vn = &.{
        .{ .fromOneBased(1), .fromOneBased(2), .fromOneBased(3) },
        .{ .fromOneBased(4), .fromOneBased(5), .fromOneBased(6) },
        .{ .fromOneBased(7), .fromOneBased(8), .fromOneBased(9) },
        .{ .fromOneBased(10), .fromOneBased(11), .fromOneBased(12) },
        .{ .fromOneBased(13), .fromOneBased(14), .fromOneBased(15) },
    } } }});

    // -- six -- //
    try testAstParse("f 1 2 3 4 5 6", &.{.{ .f = .{ .vertex_only = &.{
        .fromOneBased(1),
        .fromOneBased(2),
        .fromOneBased(3),
        .fromOneBased(4),
        .fromOneBased(5),
        .fromOneBased(6),
    } } }});
    try testAstParse("f 1/2 3/4 5/6 7/8 9/10 11/12", &.{.{ .f = .{ .only_vt = &.{
        .{ .fromOneBased(1), .fromOneBased(2) },
        .{ .fromOneBased(3), .fromOneBased(4) },
        .{ .fromOneBased(5), .fromOneBased(6) },
        .{ .fromOneBased(7), .fromOneBased(8) },
        .{ .fromOneBased(9), .fromOneBased(10) },
        .{ .fromOneBased(11), .fromOneBased(12) },
    } } }});
    try testAstParse("f 1//2 3//4 5//6 7//8 9//10 11//12", &.{.{ .f = .{ .only_vn = &.{
        .{ .fromOneBased(1), .fromOneBased(2) },
        .{ .fromOneBased(3), .fromOneBased(4) },
        .{ .fromOneBased(5), .fromOneBased(6) },
        .{ .fromOneBased(7), .fromOneBased(8) },
        .{ .fromOneBased(9), .fromOneBased(10) },
        .{ .fromOneBased(11), .fromOneBased(12) },
    } } }});
    try testAstParse("f 1/2/3 4/5/6 7/8/9 10/11/12 13/14/15 16/17/18", &.{.{ .f = .{ .both_vt_vn = &.{
        .{ .fromOneBased(1), .fromOneBased(2), .fromOneBased(3) },
        .{ .fromOneBased(4), .fromOneBased(5), .fromOneBased(6) },
        .{ .fromOneBased(7), .fromOneBased(8), .fromOneBased(9) },
        .{ .fromOneBased(10), .fromOneBased(11), .fromOneBased(12) },
        .{ .fromOneBased(13), .fromOneBased(14), .fromOneBased(15) },
        .{ .fromOneBased(16), .fromOneBased(17), .fromOneBased(18) },
    } } }});
}

test "error: null byte" {
    try testAstDiagnostic("\n\x00", .{ .InvalidBytes = .{ .line = 1, .src = "\x00" } });
}

test "error: invalid escape" {
    try testAstDiagnostic("\\a", .{ .InvalidEscape = .{ .line = 0, .src = "\\a" } });
}

test "error: invalid parameter (cstype)" {
    try testAstDiagnostic("cstype foo", .{ .InvalidParameter = .{ .line = 0, .src = "foo" } });
    try testAstDiagnostic("cstype rat\\\n foo", .{ .InvalidParameter = .{ .line = 1, .src = "foo" } });
}

test "error: invalid prefix" {
    try testAstDiagnostic("foo bar", .{ .InvalidPrefix = .{ .line = 0, .src = "foo" } });
}

test "error: invalid reference index triplet" {
    try testAstDiagnostic("f a", .{ .InvalidReferenceIndexTriplet = .{ .line = 0, .src = "a" } });
    try testAstDiagnostic("f 1/b", .{ .InvalidReferenceIndexTriplet = .{ .line = 0, .src = "1/b" } });
    try testAstDiagnostic("f 1/2/c", .{ .InvalidReferenceIndexTriplet = .{ .line = 0, .src = "1/2/c" } });
}

test "error: incomplete statement" {
    try testAstDiagnostic("v", .{ .IncompleteStatement = .{ .line = 0, .src = "v" } });
    try testAstDiagnostic("v 1", .{ .IncompleteStatement = .{ .line = 0, .src = "v" } });
    try testAstDiagnostic("v 1 2", .{ .IncompleteStatement = .{ .line = 0, .src = "v" } });
}

test "error: unexpected parameter" {
    try testAstDiagnostic("v 1 2 3 4 5", .{ .UnexpectedParameter = .{ .line = 0, .src = "5" } });
}

test "error: mixed face layouts" {
    try testAstDiagnostic("f 1/2/3 \\\n 1//3", .{ .MixedFaceLayouts = .{
        .first = .{ 0, .both_vt_vn, "1/2/3" },
        .bad = .{ 1, .only_vn, "1//3" },
    } });
    try testAstDiagnostic("f 1/2/3 \\\n 4/5/6 \\\n 1//3", .{ .MixedFaceLayouts = .{
        .first = .{ 0, .both_vt_vn, "1/2/3" },
        .bad = .{ 2, .only_vn, "1//3" },
    } });
}
