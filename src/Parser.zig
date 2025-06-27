const std = @import("std");
const stridx = @import("stridx.zig");
const Lexer = @import("Lexer.zig");
const Ast = @import("Ast.zig");

comptime {
    _ = stridx;
    _ = Lexer;
    _ = Ast;
}

const Parser = @This();
lexer: Lexer,
line: u32,
debug: if (Debug.enabled) Debug else struct {
    pub const init: @This() = .{};
},

pub const ParseOutput = struct {
    nodes: *std.MultiArrayList(Ast.Node),
    extra: *std.ArrayListUnmanaged(u32),
    indices_wip: *IndicesWip,
    str_set_wip: *StrSetWip,
};

pub fn parse(
    self: *Parser,
    gpa: std.mem.Allocator,
    maybe_diag: ?*ParseDiagnostic,
    output: ParseOutput,
) (ParseError || std.mem.Allocator.Error)!void {
    var dummy_diag: ParseDiagnostic = undefined;
    const diag = maybe_diag orelse &dummy_diag;
    diag.* = undefined;
    self.parseUnchecked(gpa, diag, output) catch |err| switch (err) {
        error.OutOfMemory => |e| return e,
        inline else => |e| {
            _ = &@field(diag, @errorName(e));
            return e;
        },
    };
}

pub const ParseError = error{
    InvalidBytes,
    InvalidEscape,
    InvalidParameter,
    InvalidPrefix,
    InvalidReferenceIndexTriplet,
    IncompleteStatement,
    UnexpectedParameter,
    MixedFaceLayouts,
};

pub const ParseDiagnostic = union {
    InvalidBytes: LineLoc,
    InvalidEscape: LineLoc,
    InvalidParameter: LineLoc,
    InvalidPrefix: LineLoc,
    InvalidReferenceIndexTriplet: LineLoc,
    IncompleteStatement: LineLoc,
    UnexpectedParameter: LineLoc,
    MixedFaceLayouts: MixedFaceLayoutsPayload,

    pub const LineLoc = struct {
        line: u32,
        loc: Lexer.Token.Loc,
    };

    pub const MixedFaceLayoutsPayload = struct {
        first: Face,
        bad: Face,

        pub const Face = struct {
            line: u32,
            loc: Lexer.Token.Loc,
            layout: Ast.FaceTriplet.Layout,
        };
    };

    pub fn tagged(self: ParseDiagnostic, err: ParseError) Tagged {
        return switch (err) {
            inline else => |e| @unionInit(
                Tagged,
                @errorName(e),
                @field(self, @errorName(e)),
            ),
        };
    }

    pub const Tagged = union(Tag) {
        pub const Tag = std.meta.FieldEnum(ParseError);
        InvalidBytes: LineLoc,
        InvalidEscape: LineLoc,
        InvalidParameter: LineLoc,
        InvalidPrefix: LineLoc,
        InvalidReferenceIndexTriplet: LineLoc,
        IncompleteStatement: LineLoc,
        UnexpectedParameter: LineLoc,
        MixedFaceLayouts: MixedFaceLayoutsPayload,

        pub fn fmt(self: Tagged, src: []const u8) Fmt {
            return .{
                .diag = self,
                .src = src,
            };
        }

        pub const Fmt = struct {
            diag: Tagged,
            src: []const u8,

            pub fn format(
                self: Fmt,
                comptime fmt_str: []const u8,
                fmt_options: std.fmt.FormatOptions,
                writer: anytype,
            ) @TypeOf(writer).Error!void {
                _ = fmt_str;
                _ = fmt_options;

                switch (self.diag) {
                    inline //
                    .InvalidBytes,
                    .InvalidEscape,
                    .InvalidParameter,
                    .InvalidPrefix,
                    .InvalidReferenceIndexTriplet,
                    .UnexpectedParameter,
                    .IncompleteStatement,
                    => |line_loc, tag| {
                        const str = line_loc.loc.getStr(self.src);
                        try writer.print("{s} on line {d}: '{s}'", .{ @tagName(tag), line_loc.line, str });
                    },

                    inline .MixedFaceLayouts => |mfl, e| {
                        try writer.print("{s}: '{s}' on line {d} vs '{s}' on line {d}", .{
                            @tagName(e),
                            mfl.first.loc.getStr(self.src),
                            mfl.first.line,
                            mfl.bad.loc.getStr(self.src),
                            mfl.bad.line,
                        });
                    },
                }
            }
        };
    };

    comptime {
        const e_info = @typeInfo(ParseError).error_set.?;
        const u_info = @typeInfo(ParseDiagnostic).@"union";
        const amt = @min(e_info.len, u_info.fields.len);
        const idx_fields = @typeInfo(@TypeOf(.{{}} ** amt)).@"struct".fields;
        for (e_info[0..amt], u_info.fields[0..amt], idx_fields) |e_field, u_field, idx_field| {
            const mismatch_err = "[" ++ idx_field.name ++ "] " ++ "Expected `" ++ e_field.name ++ "`, instead found `" ++ u_field.name ++ "`";
            if (!ctEql(e_field.name, u_field.name)) @compileError(mismatch_err);
        }
        switch (std.math.order(e_info.len, u_info.fields.len)) {
            .eq => {},
            .lt => {
                var msg: []const u8 = "Error missing fields: ";
                for (u_info.fields[amt..], 0..) |u_field, i| {
                    if (i != 0) msg = msg ++ ", ";
                    msg = msg ++ u_field.name;
                }
                @compileError(msg);
            },
            .gt => {
                var msg: []const u8 = "Union missing fields: ";
                for (e_info[amt..], 0..) |e_field, i| {
                    if (i != 0) msg = msg ++ ", ";
                    msg = msg ++ e_field.name;
                }
                @compileError(msg);
            },
        }
    }

    inline fn ctEql(comptime a: []const u8, comptime b: []const u8) bool {
        comptime {
            if (a.len != b.len) return false;
            const e_name: @Vector(a.len, u8) = a[0..].*;
            const u_name: @Vector(b.len, u8) = b[0..].*;
            return @reduce(.And, e_name == u_name);
        }
    }

    inline fn ret(
        self: *ParseDiagnostic,
        comptime err: ParseError,
        data: @FieldType(ParseDiagnostic, @errorName(err)),
    ) @TypeOf(switch (err) {
        err => |e| e,
        else => unreachable,
    }) {
        self.* = @unionInit(ParseDiagnostic, @errorName(err), data);
        return switch (err) {
            err => |e| e, // weird hack for narrowing the comptime-known error's type
            else => comptime unreachable,
        };
    }
};

fn parseUnchecked(
    self: *Parser,
    gpa: std.mem.Allocator,
    diag: *ParseDiagnostic,
    output: ParseOutput,
) (ParseError || std.mem.Allocator.Error)!void {
    const src = self.lexer.src;
    std.debug.assert(src.len <= std.math.maxInt(u32));

    const nodes = output.nodes;
    const extra = output.extra;
    const indices_wip = output.indices_wip;
    const str_set_wip = output.str_set_wip;

    mainloop: while (true) {
        const prefix_line = self.line;
        const maybe_prefix: ?Ast.Node.Full.Prefix, const prefix_loc: Lexer.Token.Loc = switch (self.nextTok()) {
            .str => |first_loc| blk: {
                var str_iter = self.stringIter(first_loc) orelse break :blk .{
                    std.meta.stringToEnum(Ast.Node.Full.Prefix, first_loc.getStr(src)),
                    first_loc,
                };
                const maybe_prefix = str_iter.peekToEnum(Ast.Node.Full.Prefix);
                str_iter.consume();
                break :blk .{ maybe_prefix, str_iter.fullSrcLoc() };
            },
            .invalid => |loc| return diag.ret(error.InvalidBytes, .{
                .line = self.line,
                .loc = loc,
            }),
            .backslash_invalid => |loc| return diag.ret(error.InvalidEscape, .{
                .line = self.line,
                .loc = .{
                    .start = loc.start,
                    .end = loc.end + @intFromBool(loc.end != self.lexer.src.len),
                },
            }),

            .nl => continue :mainloop,
            .eof => break :mainloop,

            .comment_nl_one => continue :mainloop,
            .comment_nl_two => continue :mainloop,
            .comment_eof => continue :mainloop,

            .backslash_nl => continue :mainloop,
            .whitespace => continue :mainloop,
        };

        const prefix = maybe_prefix orelse {
            return diag.ret(error.InvalidPrefix, .{
                .line = prefix_line,
                .loc = prefix_loc,
            });
        };

        switch (prefix) {
            .v => {
                try indices_wip.ensureUnusedCapacity(gpa, 1);
                try nodes.ensureUnusedCapacity(gpa, 1);
                try extra.ensureUnusedCapacity(gpa, 4);
                const node_index: Ast.Node.Index = @intCast(nodes.len);
                const extra_index: u32 = @intCast(extra.items.len);

                const v_x = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                    .str_idx => |str_idx| str_idx,
                    .line_end => return diag.ret(error.IncompleteStatement, .{
                        .line = prefix_line,
                        .loc = prefix_loc,
                    }),
                };
                const v_y = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                    .str_idx => |str_idx| str_idx,
                    .line_end => return diag.ret(error.IncompleteStatement, .{
                        .line = prefix_line,
                        .loc = prefix_loc,
                    }),
                };
                const v_z = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                    .str_idx => |str_idx| str_idx,
                    .line_end => return diag.ret(error.IncompleteStatement, .{
                        .line = prefix_line,
                        .loc = prefix_loc,
                    }),
                };
                const v_w = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                    .str_idx => |str_idx| str_idx,
                    .line_end => {
                        const data: *[3]stridx.Index = elemCast(stridx.Index, extra.addManyAsArrayAssumeCapacity(3));
                        data.* = .{ v_x, v_y, v_z };
                        nodes.appendAssumeCapacity(.{ .tag = .v_xyz, .data = .{ .index = extra_index } });
                        indices_wip.addAssumeCapacity(.positions, node_index);
                        continue :mainloop;
                    },
                };
                try self.expectEndOfParameterList(diag);

                const data: *[4]stridx.Index = elemCast(stridx.Index, extra.addManyAsArrayAssumeCapacity(4));
                data.* = .{ v_x, v_y, v_z, v_w };
                nodes.appendAssumeCapacity(.{ .tag = .v_xyzw, .data = .{ .index = extra_index } });
                indices_wip.addAssumeCapacity(.positions, node_index);
                continue :mainloop;
            },

            .vn => {
                try indices_wip.ensureUnusedCapacity(gpa, 1);
                try nodes.ensureUnusedCapacity(gpa, 1);
                try extra.ensureUnusedCapacity(gpa, 3);
                const node_index: Ast.Node.Index = @intCast(nodes.len);
                const extra_index: u32 = @intCast(extra.items.len);

                const vn_i = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                    .str_idx => |str_idx| str_idx,
                    .line_end => return diag.ret(error.IncompleteStatement, .{
                        .line = prefix_line,
                        .loc = prefix_loc,
                    }),
                };
                const vn_j = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                    .str_idx => |str_idx| str_idx,
                    .line_end => return diag.ret(error.IncompleteStatement, .{
                        .line = prefix_line,
                        .loc = prefix_loc,
                    }),
                };
                const vn_k = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                    .str_idx => |str_idx| str_idx,
                    .line_end => return diag.ret(error.IncompleteStatement, .{
                        .line = prefix_line,
                        .loc = prefix_loc,
                    }),
                };
                try self.expectEndOfParameterList(diag);

                const data: *[3]stridx.Index = elemCast(stridx.Index, extra.addManyAsArrayAssumeCapacity(3));
                data.* = .{ vn_i, vn_j, vn_k };
                nodes.appendAssumeCapacity(.{ .tag = .vn_ijk, .data = .{ .index = extra_index } });
                indices_wip.addAssumeCapacity(.normals, node_index);
                continue :mainloop;
            },

            inline .vt, .vp => |tag| {
                const u: Ast.Node.Tag, //
                const uv: Ast.Node.Tag, //
                const uvw: Ast.Node.Tag, //
                const indices_offset_field: Ast.Indices.Offsets.FieldTag //
                = comptime switch (tag) {
                    .vt => .{ .vt_u, .vt_uv, .vt_uvw, .textures },
                    .vp => .{ .vp_u, .vp_uv, .vp_uvw, .parameter_spaces },
                    else => unreachable,
                };

                try indices_wip.ensureUnusedCapacity(gpa, 1);
                try nodes.ensureUnusedCapacity(gpa, 1);
                try extra.ensureUnusedCapacity(gpa, 3);
                const node_index: Ast.Node.Index = @intCast(nodes.len);
                const extra_index: u32 = @intCast(extra.items.len);

                const vp_u = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                    .str_idx => |str_idx| str_idx,
                    .line_end => return diag.ret(error.IncompleteStatement, .{
                        .line = prefix_line,
                        .loc = prefix_loc,
                    }),
                };
                const vp_uv = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                    .str_idx => |str_idx| str_idx,
                    .line_end => {
                        nodes.appendAssumeCapacity(.{ .tag = u, .data = .{ .value = vp_u } });
                        indices_wip.addAssumeCapacity(indices_offset_field, node_index);
                        continue :mainloop;
                    },
                };
                const vp_uvw = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                    .str_idx => |str_idx| str_idx,
                    .line_end => {
                        const data: *[2]stridx.Index = elemCast(stridx.Index, extra.addManyAsArrayAssumeCapacity(2));
                        data.* = .{ vp_u, vp_uv };
                        nodes.appendAssumeCapacity(.{ .tag = uv, .data = .{ .index = extra_index } });
                        indices_wip.addAssumeCapacity(indices_offset_field, node_index);
                        continue :mainloop;
                    },
                };
                try self.expectEndOfParameterList(diag);

                const data: *[3]stridx.Index = elemCast(stridx.Index, extra.addManyAsArrayAssumeCapacity(3));
                data.* = .{ vp_u, vp_uv, vp_uvw };
                nodes.appendAssumeCapacity(.{ .tag = uvw, .data = .{ .index = extra_index } });
                indices_wip.addAssumeCapacity(indices_offset_field, node_index);
                continue :mainloop;
            },

            .cstype => {
                try nodes.ensureUnusedCapacity(gpa, 1);

                switch (try self.expectSpaceSepAfterString(diag)) {
                    .line_end => return diag.ret(error.IncompleteStatement, .{
                        .line = prefix_line,
                        .loc = prefix_loc,
                    }),
                    .whitespace => {},
                }
                const rat_first_loc = switch (try self.expectStringFirst(diag)) {
                    .first_loc => |first_loc| first_loc,
                    .line_end => return diag.ret(error.IncompleteStatement, .{
                        .line = prefix_line,
                        .loc = prefix_loc,
                    }),
                    .whitespace => unreachable,
                };
                var maybe_rat_str_iter = self.stringIter(rat_first_loc);
                const is_rational = if (tryStrToEnum(
                    self.lexer.src,
                    rat_first_loc,
                    maybe_rat_str_iter,
                    enum { rat },
                )) |rat| switch (rat) {
                    .rat => blk: {
                        if (maybe_rat_str_iter) |*str_iter| str_iter.consume();
                        break :blk true;
                    },
                } else false;

                if (is_rational) switch (try self.expectSpaceSepAfterString(diag)) {
                    .line_end => return diag.ret(error.IncompleteStatement, .{
                        .line = prefix_line,
                        .loc = prefix_loc,
                    }),
                    .whitespace => {},
                };

                const cstype_tag: Ast.Node.CsType.Type = cstype_tag: {
                    const type_loc, var maybe_type_str_iter = if (!is_rational)
                        .{ rat_first_loc, maybe_rat_str_iter }
                    else switch (try self.expectStringFirst(diag)) {
                        .first_loc => |type_first_loc| .{ type_first_loc, self.stringIter(type_first_loc) },
                        .line_end => return diag.ret(error.IncompleteStatement, .{
                            .line = prefix_line,
                            .loc = prefix_loc,
                        }),
                        .whitespace => unreachable,
                    };

                    const cstype_tag_line = self.line;
                    const cstype_tag = tryStrToEnum(
                        self.lexer.src,
                        type_loc,
                        maybe_type_str_iter,
                        Ast.Node.CsType.Type,
                    ) orelse return diag.ret(error.InvalidParameter, .{
                        .line = cstype_tag_line,
                        .loc = type_loc,
                    });
                    if (maybe_type_str_iter) |*str_iter| str_iter.consume();
                    break :cstype_tag cstype_tag;
                };

                const cstype: Ast.Node.CsType = .{
                    .rational = is_rational,
                    .type = cstype_tag,
                };
                nodes.appendAssumeCapacity(.{ .tag = .cstype, .data = .{ .cstype = .{ .value = cstype } } });
                continue :mainloop;
            },

            .deg => {
                try nodes.ensureUnusedCapacity(gpa, 1);
                try extra.ensureUnusedCapacity(gpa, 2);
                const extra_index: u32 = @intCast(extra.items.len);

                const degu_str_idx = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                    .str_idx => |str_idx| str_idx,
                    .line_end => return diag.ret(error.IncompleteStatement, .{
                        .line = prefix_line,
                        .loc = prefix_loc,
                    }),
                };

                const degv_str_idx_opt = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                    .str_idx => |str_idx| blk: {
                        try self.expectEndOfParameterList(diag);
                        break :blk str_idx;
                    },
                    .line_end => null,
                };

                if (degv_str_idx_opt) |degv_str_idx| {
                    const data: *[2]stridx.Index = elemCast(stridx.Index, extra.addManyAsArrayAssumeCapacity(2));
                    data.* = .{ degu_str_idx, degv_str_idx };
                    nodes.appendAssumeCapacity(.{ .tag = .deg_surface, .data = .{ .index = extra_index } });
                } else {
                    nodes.appendAssumeCapacity(.{ .tag = .deg_curve, .data = .{ .value = degu_str_idx } });
                }

                continue :mainloop;
            },

            inline .o, .s, .usemtl => |tag| {
                try nodes.ensureUnusedCapacity(gpa, 1);

                const str_idx = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                    .str_idx => |str_idx| str_idx,
                    .line_end => return diag.ret(error.IncompleteStatement, .{
                        .line = prefix_line,
                        .loc = prefix_loc,
                    }),
                };

                nodes.appendAssumeCapacity(.{
                    .tag = @field(Ast.Node.Tag, @tagName(tag)),
                    .data = .{ .string = str_idx },
                });

                continue :mainloop;
            },

            inline .g, .call, .mtllib => |tag| {
                try nodes.ensureUnusedCapacity(gpa, 1);
                const extra_index: u32 = @intCast(extra.items.len);

                const node_tag_empty: ?Ast.Node.Tag, //
                const node_tag_single: Ast.Node.Tag, //
                const node_tag_multi: Ast.Node.Tag //
                = comptime switch (tag) {
                    // zig fmt: off
                    .g =>      .{ .g_empty, .g_single,   .g_multi },
                    .call =>   .{ null,     .call,       .call_no_args },
                    .mtllib => .{ null,     .mtllib_one, .mtllib_multi },
                    else => unreachable,
                    // zig fmt: on
                };

                const first_str_idx = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                    .str_idx => |str_idx| str_idx,
                    .line_end => {
                        const node_tag = node_tag_empty orelse return diag.ret(error.IncompleteStatement, .{
                            .line = prefix_line,
                            .loc = prefix_loc,
                        });
                        nodes.appendAssumeCapacity(.{ .tag = node_tag, .data = .{ .empty = {} } });
                        continue :mainloop;
                    },
                };

                while (true) {
                    const str_idx = switch (try self.expectSpaceAndParamStrAfterString(gpa, diag, str_set_wip)) {
                        .str_idx => |str_idx| str_idx,
                        .line_end => break,
                    };

                    if (extra_index == extra.items.len) {
                        try extra.ensureUnusedCapacity(gpa, 2);
                        extra.appendAssumeCapacity(first_str_idx.value().?);
                    }

                    try extra.ensureUnusedCapacity(gpa, 2);
                    extra.appendAssumeCapacity(str_idx.value().?);
                }

                if (extra_index == extra.items.len) {
                    nodes.appendAssumeCapacity(.{ .tag = node_tag_single, .data = .{ .value = first_str_idx } });
                } else {
                    extra.appendAssumeCapacity(stridx.Index.valueAllowNull(.null));
                    nodes.appendAssumeCapacity(.{ .tag = node_tag_multi, .data = .{ .index = extra_index } });
                }

                continue :mainloop;
            },

            .f => {
                try nodes.ensureUnusedCapacity(gpa, 1);
                const node_index: u32 = @intCast(nodes.len);
                const extra_index: u32 = @intCast(extra.items.len);

                switch (try self.expectSpaceSepAfterString(diag)) {
                    .line_end => return diag.ret(error.IncompleteStatement, .{
                        .line = prefix_line,
                        .loc = prefix_loc,
                    }),
                    .whitespace => {},
                }

                const first_triplet_line = self.line;
                const first_triplet = try self.expectRefIdxTriplet(indices_wip.asIndices(), diag) orelse {
                    return diag.ret(error.IncompleteStatement, .{
                        .line = prefix_line,
                        .loc = prefix_loc,
                    });
                };
                const first_tok_end: Parser.TokEnd =
                    switch (try self.expectSpaceSepAfterString(diag)) {
                        .line_end => .line_end,
                        .whitespace => .whitespace,
                    };
                const first_layout = first_triplet.value.layout();
                try zipUpFaceTriplet(gpa, extra, first_triplet.value, switch (first_tok_end) {
                    .line_end => .no_reserve_capacity_for_null_sentinel,
                    .whitespace => .do_reserve_capacity_for_null_sentinel,
                });

                const node_tag_one: Ast.Node.Tag, //
                const node_tag_multi: Ast.Node.Tag //
                = switch (first_layout) {
                    // zig fmt: off
                    .vertex_only => .{ .f_one_v,       .f_multi_v       },
                    .only_vt     => .{ .f_one_v_vt,    .f_multi_v_vt    },
                    .only_vn     => .{ .f_one_v_vn,    .f_multi_v_vn    },
                    .both_vt_vn  => .{ .f_one_v_vt_vn, .f_multi_v_vt_vn },
                    // zig fmt: on
                };

                switch (first_tok_end) {
                    .line_end => {
                        nodes.appendAssumeCapacity(.{
                            .tag = node_tag_one,
                            .data = if (first_layout == .vertex_only)
                                .{ .ref_index = first_triplet.value.v.nonNull() }
                            else
                                .{ .index = extra_index },
                        });
                        continue :mainloop;
                    },
                    .whitespace => {
                        nodes.appendAssumeCapacity(.{
                            .tag = node_tag_multi,
                            .data = .{ .index = extra_index },
                        });
                    },
                }

                var f_triplet_count: u32 = 1;
                while (true) {
                    defer f_triplet_count += 1;
                    self.debug.assertPrevTokenKindWas(.whitespace);

                    const curr_triplet_line = self.line;
                    const curr_triplet = try self.expectRefIdxTriplet(indices_wip.asIndices(), diag) orelse {
                        return diag.ret(error.IncompleteStatement, .{
                            .line = prefix_line,
                            .loc = prefix_loc,
                        });
                    };
                    const curr_tok_end: Parser.TokEnd =
                        switch (try self.expectSpaceSepAfterString(diag)) {
                            .line_end => .line_end,
                            .whitespace => .whitespace,
                        };
                    const curr_layout = curr_triplet.value.layout();

                    if (curr_layout != first_layout) {
                        return diag.ret(error.MixedFaceLayouts, .{
                            // zig fmt: off
                            .first = .{ .line = first_triplet_line, .loc = first_triplet.loc, .layout = first_layout },
                            .bad   = .{ .line =  curr_triplet_line, .loc =  curr_triplet.loc, .layout =  curr_layout },
                            // zig fmt: on
                        });
                    }

                    try zipUpFaceTriplet(gpa, extra, curr_triplet.value, .do_reserve_capacity_for_null_sentinel);

                    switch (curr_tok_end) {
                        .line_end => break,
                        .whitespace => {},
                    }
                }

                const tag: *Ast.Node.Tag = &nodes.items(.tag)[node_index];
                switch (f_triplet_count) {
                    0, 1 => unreachable,
                    2 => tag.* = switch (tag.*) {
                        .f_multi_v => .f_two_v,
                        .f_multi_v_vt => .f_two_v_vt,
                        .f_multi_v_vn => .f_two_v_vn,
                        .f_multi_v_vt_vn => .f_two_v_vt_vn,
                        else => unreachable,
                    },
                    3 => tag.* = switch (tag.*) {
                        .f_multi_v => .f_three_v,
                        .f_multi_v_vt => .f_three_v_vt,
                        .f_multi_v_vn => .f_three_v_vn,
                        .f_multi_v_vt_vn => .f_three_v_vt_vn,
                        else => unreachable,
                    },
                    4 => tag.* = switch (tag.*) {
                        .f_multi_v => .f_four_v,
                        .f_multi_v_vt => .f_four_v_vt,
                        .f_multi_v_vn => .f_four_v_vn,
                        .f_multi_v_vt_vn => .f_four_v_vt_vn,
                        else => unreachable,
                    },
                    else => extra.appendAssumeCapacity(
                        @bitCast(Ast.RefIndex.valueAllowNull(.null)),
                    ),
                }

                continue :mainloop;
            },
        }
    }
}

pub const IndicesWip = struct {
    list: std.ArrayListUnmanaged(Ast.Node.Index),
    offsets: Ast.Indices.Offsets,

    pub const empty: IndicesWip = .{
        .list = .empty,
        .offsets = .empty,
    };

    pub fn deinit(self: IndicesWip, gpa: std.mem.Allocator) void {
        var list = self.list;
        list.deinit(gpa);
    }

    /// Transfers ownership of the list to the caller.
    /// Makes `self.deinit(gpa)` a noop.
    pub fn toOwnedIndices(self: *IndicesWip, gpa: std.mem.Allocator) std.mem.Allocator.Error!Ast.Indices {
        return .{
            .list = try self.list.toOwnedSlice(gpa),
            .offsets = self.offsets,
        };
    }

    /// Returns a view of the WIP indices. It is illegal to deinit the returned view.
    fn asIndices(self: IndicesWip) Ast.Indices {
        return .{
            .list = self.list.items,
            .offsets = self.offsets,
        };
    }

    fn ensureUnusedCapacity(
        self: *IndicesWip,
        gpa: std.mem.Allocator,
        additional_count: usize,
    ) std.mem.Allocator.Error!void {
        try self.list.ensureUnusedCapacity(gpa, additional_count);
    }

    fn addAssumeCapacity(
        self: *IndicesWip,
        comptime field: Ast.Indices.Offsets.FieldTag,
        node_index: Ast.Node.Index,
    ) void {
        const insert_pos: u32 = pos: {
            const end_tag = (comptime Ast.Indices.Offsets.endFieldOf(field)) orelse
                break :pos @intCast(self.list.items.len);
            break :pos @field(self.offsets, @tagName(end_tag));
        };
        self.list.insertAssumeCapacity(insert_pos, node_index);

        comptime var start_tag = field;
        inline while (true) {
            const end_tag = comptime Ast.Indices.Offsets.endFieldOf(start_tag) orelse break;
            defer start_tag = end_tag;
            @field(self.offsets, @tagName(end_tag)) += 1;
        }
    }
};

pub const StrSetWip = struct {
    bytes: std.ArrayListUnmanaged(u8),
    indexer: stridx.Indexer,

    pub const empty: StrSetWip = .{
        .bytes = .empty,
        .indexer = .empty,
    };

    pub fn deinit(self: StrSetWip, gpa: std.mem.Allocator) void {
        var bytes = self.bytes;
        bytes.deinit(gpa);
        var indexer = self.indexer;
        indexer.deinit(gpa);
    }

    /// Transfers ownership of the string set to the caller.
    /// Makes `self.deinit(gpa)` a noop.
    pub fn toOwnedStrSet(self: *StrSetWip, gpa: std.mem.Allocator) std.mem.Allocator.Error!Ast.StrSet {
        return .{
            .bytes = try self.bytes.toOwnedSlice(gpa),
            .indexer = self.indexer.move(),
        };
    }

    fn getOrPutStr(
        self: *StrSetWip,
        gpa: std.mem.Allocator,
        str: []const u8,
    ) std.mem.Allocator.Error!stridx.Index {
        std.debug.assert(std.mem.indexOfScalar(u8, str, 0) == null);

        try self.indexer.ensureUnusedCapacityContext(gpa, 1, stridx.hashCtx(self.bytes.items));
        try self.bytes.ensureUnusedCapacity(gpa, str.len + 1);
        const gop = self.indexer.getOrPutAssumeCapacityAdapted(
            str,
            stridx.hashStrAdapter(self.bytes.items),
        );
        if (!gop.found_existing) {
            gop.key_ptr.* = .from(@intCast(self.bytes.items.len));
            gop.value_ptr.* = {};
            self.bytes.appendSliceAssumeCapacity(str);
            self.bytes.appendAssumeCapacity(0);
        }
        return gop.key_ptr.*;
    }

    fn getOrPutStringIter(
        self: *StrSetWip,
        gpa: std.mem.Allocator,
        str_iter: StringIter,
    ) std.mem.Allocator.Error!stridx.Index {
        std.debug.assert(str_iter.state == .first_pending);

        try self.indexer.ensureUnusedCapacityContext(gpa, 1, stridx.hashCtx(self.bytes.items));
        try self.bytes.ensureUnusedCapacity(gpa, str_iter.peekLength() + 1);

        const ResettableIter = struct {
            iter: *const StringIter,
            current_peeker: *StringIter,

            pub fn reset(iter: @This()) void {
                iter.current_peeker.* = iter.iter.toPeeker();
            }

            pub fn next(iter: @This()) ?[]const u8 {
                const segment, _ = iter.current_peeker.next() orelse return null;
                return segment;
            }
        };
        var current_peeker = str_iter.toPeeker();
        const gop = self.indexer.getOrPutAssumeCapacityAdapted(
            ResettableIter{ .iter = &str_iter, .current_peeker = &current_peeker },
            stridx.hashIterAdapter(self.bytes.items, ResettableIter),
        );
        if (!gop.found_existing) {
            gop.key_ptr.* = .from(@intCast(self.bytes.items.len));
            gop.value_ptr.* = {};
            var peeker = str_iter.toPeeker();
            while (true) {
                const segment, _ = peeker.next() orelse break;
                self.bytes.appendSliceAssumeCapacity(segment);
            }
            self.bytes.appendAssumeCapacity(0);
        }
        return gop.key_ptr.*;
    }
};

const Debug = struct {
    const enabled = @import("builtin").mode == .Debug;
    state: if (Debug.enabled) RealState else EmptyState,

    pub const init: Debug = .{ .state = .init };

    const RealState = struct {
        prev_tok: ?Lexer.Token,

        const init: RealState = .{
            .prev_tok = null,
        };
    };

    const EmptyState = struct {
        const init: EmptyState = .{};
    };

    fn assertPrevTokenEql(self: *const Debug, expected: Lexer.Token) void {
        if (!Debug.enabled) return;
        const actual = self.state.prev_tok orelse std.debug.panic("Expected {}, found null", .{expected});
        if (!expected.eql(actual)) {
            std.debug.panic("Expected {}, found {}", .{ expected, actual });
        }
    }
    fn assertPrevTokenKindWasNot(
        self: *const Debug,
        sentinel_kind: Lexer.Token.Kind,
    ) void {
        if (!Debug.enabled) return;
        const actual = self.state.prev_tok orelse return;
        if (actual == sentinel_kind) std.debug.panic(
            "Expected something other than {s}",
            .{@tagName(sentinel_kind)},
        );
    }
    fn assertPrevTokenKindWas(
        self: *const Debug,
        expected_kind: Lexer.Token.Kind,
    ) void {
        if (!Debug.enabled) return;
        const actual = self.state.prev_tok orelse std.debug.panic("Expected {s}, found null", .{@tagName(expected_kind)});
        if (actual != expected_kind) std.debug.panic(
            "Expected {s}, found {s}",
            .{ @tagName(expected_kind), @tagName(actual) },
        );
    }
    fn assertPrevTokenKindWasOneOf(
        self: *const Debug,
        kinds: []const Lexer.Token.Kind,
    ) void {
        if (!Debug.enabled) return;
        const actual = self.state.prev_tok orelse std.debug.panic("Expected one of {any}, found null", .{kinds});
        for (kinds) |kind| {
            if (actual == kind) break;
        } else {
            std.debug.panic("Expected one of {any}, found {}", .{ kinds, actual });
        }
    }
};

fn nextTok(self: *Parser) Lexer.Token {
    const tok = self.lexer.next();
    switch (tok) {
        .backslash_nl,
        .comment_nl_one,
        .comment_nl_two,
        .nl,
        => self.line += 1,
        else => {},
    }
    if (Debug.enabled) {
        self.debug.state.prev_tok = tok;
    }
    return tok;
}

/// Asserts the previous token was a `.backslash_nl`.
/// Checks if the next token is also an escaped newline, and if it is,
/// skips over it, and checks the next token, returning the first token
/// which isn't an escaped newline.
/// Return value will never be `.backslash_nl`.
fn nextSkipEscapedNewlineChain(self: *Parser) Lexer.Token {
    self.debug.assertPrevTokenKindWas(.backslash_nl);
    return sw: switch (self.nextTok()) {
        .backslash_nl => continue :sw self.nextTok(),
        else => |tok| tok,
    };
}

const ExpectEndOfLineError = error{
    InvalidBytes,
    InvalidEscape,
};

/// Returns null on end of line.
/// Returns the location of the unexpected token if encountered.
fn expectEndOfLine(
    self: *Parser,
    diag: *ParseDiagnostic,
) ExpectEndOfLineError!?Lexer.Token.Loc {
    sw: switch (self.nextTok()) {
        .invalid => |loc| return diag.ret(error.InvalidBytes, .{
            .line = self.line,
            .loc = loc,
        }),
        .backslash_invalid => |loc| return diag.ret(error.InvalidEscape, .{
            .line = self.line,
            .loc = .{
                .start = loc.start,
                .end = loc.end + @intFromBool(loc.end != self.lexer.src.len),
            },
        }),

        .str => |loc| return loc,
        .backslash_nl => continue :sw self.nextSkipEscapedNewlineChain(),
        .whitespace => {
            const tok = self.nextTok();
            std.debug.assert(tok != .whitespace);
            continue :sw tok;
        },

        .comment_nl_one,
        .comment_nl_two,
        .comment_eof,
        .nl,
        .eof,
        => return null,
    }
}

const ExpectEndOfParameterListError =
    ExpectEndOfLineError ||
    error{
        UnexpectedParameter,
    };

/// Expects the end of a parameter list with definite length, returning
/// `UnexpectedParameter` if there's anything other than a line end.
fn expectEndOfParameterList(
    self: *Parser,
    diag: *ParseDiagnostic,
) ExpectEndOfParameterListError!void {
    if (try self.expectEndOfLine(diag)) |unexpected_loc| {
        return diag.ret(error.UnexpectedParameter, .{
            .line = self.line,
            .loc = unexpected_loc,
        });
    }
}

const ExpectSpaceSepAfterStringError = error{
    InvalidBytes,
    InvalidEscape,
};

/// Used to expect whitespace immediately after a string, like a statement prefix.
/// If there is a newline immediately following a space, it is eagerly consumed, and returned (`.line_end`).
/// This means that this only returns `whitespace` if the whitespace is a separator for a subsequent token.
/// Escaped newlines are automatically collapsed and ignored in this process.
fn expectSpaceSepAfterString(
    self: *Parser,
    diag: *ParseDiagnostic,
) ExpectSpaceSepAfterStringError!union(enum) {
    whitespace,
    line_end,
} {
    self.debug.assertPrevTokenKindWas(.str);
    sw: switch (self.nextTok()) {
        .str => unreachable, // would have been part of the previous string

        .invalid => |loc| return diag.ret(error.InvalidBytes, .{
            .line = self.line,
            .loc = loc,
        }),
        .backslash_invalid => |loc| return diag.ret(error.InvalidEscape, .{
            .line = self.line,
            .loc = .{
                .start = loc.start,
                .end = loc.end + @intFromBool(loc.end != self.lexer.src.len),
            },
        }),

        .comment_nl_one,
        .comment_nl_two,
        .comment_eof,
        .nl,
        .eof,
        => return .line_end,

        .backslash_nl => continue :sw self.nextSkipEscapedNewlineChain(),

        .whitespace => {
            var peeker = self.*;
            peeker.debug.assertPrevTokenKindWas(.whitespace);
            peek: switch (peeker.nextTok().getKind()) {
                // we just saw a whitespace, and the previous iteration would
                // have broke the loop if it was a whitespace or a newline.
                .whitespace => unreachable,

                // we need to skip the escaped newline, and any subsequent
                // escaped newlines if there's a chain of them.
                .backslash_nl => switch (peeker.nextSkipEscapedNewlineChain().getKind()) {
                    // we just skipped any potential immediate escaped newlines
                    .backslash_nl => unreachable,
                    .whitespace => {
                        std.debug.assert(self.nextTok() == .backslash_nl);
                        std.debug.assert(self.nextSkipEscapedNewlineChain() == .whitespace);
                        // we're back to where we started: we just encountered a whitespace,
                        // and we want to peek at what's ahead to see if there's a line_end
                        // we must eagerly consume and report.
                        continue :peek peeker.nextTok();
                    },
                    else => |tag| continue :peek tag,
                },

                // if it's a newline or eof, return line_end, consuming the token.
                .comment_nl_one,
                .comment_nl_two,
                .nl,
                .comment_eof,
                .eof,
                => |tag| {
                    std.debug.assert(self.nextTok() == tag);
                    return .line_end;
                },

                // if it's anything else, just return the latest whitespace
                // without consuming the token we just peeked.
                .backslash_invalid,
                .invalid,
                .str,
                => return .whitespace,
            }
        },
    }
}

const ExpectStringFirstError = error{
    InvalidBytes,
    InvalidEscape,
};

fn expectStringFirst(
    self: *Parser,
    diag: *ParseDiagnostic,
) ExpectStringFirstError!union(enum) {
    first_loc: Lexer.Token.Loc,
    whitespace: Lexer.Token.Loc,
    line_end,
} {
    self.debug.assertPrevTokenKindWasNot(.str);
    sw: switch (self.nextTok()) {
        .invalid => |loc| return diag.ret(error.InvalidBytes, .{
            .line = self.line,
            .loc = loc,
        }),
        .backslash_invalid => |loc| return diag.ret(error.InvalidEscape, .{
            .line = self.line,
            .loc = .{
                .start = loc.start,
                .end = loc.end + @intFromBool(loc.end != self.lexer.src.len),
            },
        }),
        .whitespace => |loc| return .{ .whitespace = loc },

        .comment_eof,
        .comment_nl_one,
        .comment_nl_two,
        .nl,
        .eof,
        => return .line_end,

        .backslash_nl => continue :sw self.nextTok(),
        .str => |first_loc| return .{ .first_loc = first_loc },
    }
}

const ExpectParamStrError =
    ExpectSpaceSepAfterStringError ||
    std.mem.Allocator.Error;

fn expectParamStr(
    self: *Parser,
    gpa: std.mem.Allocator,
    diag: *ParseDiagnostic,
    str_set_wip: *StrSetWip,
) ExpectParamStrError!union(enum) {
    str_idx: stridx.Index,
    whitespace: Lexer.Token.Loc,
    line_end,
} {
    self.debug.assertPrevTokenKindWasNot(.str);
    const first_loc = switch (try self.expectStringFirst(diag)) {
        .first_loc => |first_loc| first_loc,
        .whitespace => |loc| return .{ .whitespace = loc },
        .line_end => return .line_end,
    };
    var str_iter = self.stringIter(first_loc) orelse {
        const str_idx = try str_set_wip.getOrPutStr(gpa, first_loc.getStr(self.lexer.src));
        return .{ .str_idx = str_idx };
    };
    const str_idx = try str_set_wip.getOrPutStringIter(gpa, str_iter);
    str_iter.consume();
    return .{ .str_idx = str_idx };
}

const ExpectSpaceAndParamStrAfterStringError =
    ExpectParamStrError;

/// Expects a whitespace, and then a parameter string, which
/// is then interned, and has its interned string index returned.
fn expectSpaceAndParamStrAfterString(
    self: *Parser,
    gpa: std.mem.Allocator,
    diag: *ParseDiagnostic,
    str_set_wip: *StrSetWip,
) ExpectSpaceAndParamStrAfterStringError!union(enum) {
    str_idx: stridx.Index,
    line_end,
} {
    switch (try self.expectSpaceSepAfterString(diag)) {
        .line_end => return .line_end,
        .whitespace => {},
    }

    switch (try self.expectParamStr(gpa, diag, str_set_wip)) {
        .str_idx => |str_idx| return .{ .str_idx = str_idx },
        .whitespace => unreachable,
        .line_end => unreachable,
    }
}

const TokEnd = enum {
    whitespace,
    line_end,
};

const RefIdxTripletResult = struct {
    loc: Lexer.Token.Loc,
    value: Ast.FaceTriplet,
};

const ExpectRefIdxTripletError =
    ExpectStringFirstError ||
    error{
        InvalidReferenceIndexTriplet,
    };

/// Returns null on immediate line end.
fn expectRefIdxTriplet(
    self: *Parser,
    indices: Ast.Indices,
    diag: *ParseDiagnostic,
) ExpectRefIdxTripletError!?RefIdxTripletResult {
    self.debug.assertPrevTokenKindWas(.whitespace);

    const first_loc, const line: u32 = sw: switch (self.nextTok()) {
        .str => |first_loc| .{ first_loc, self.line },
        .backslash_nl => continue :sw self.nextSkipEscapedNewlineChain(),
        // we just saw a whitespace
        .whitespace => unreachable,

        .comment_nl_one,
        .comment_nl_two,
        .nl,
        .comment_eof,
        .eof,
        => return null,

        .backslash_invalid => |loc| return diag.ret(error.InvalidEscape, .{
            .line = self.line,
            .loc = .{
                .start = loc.start,
                .end = loc.end + @intFromBool(loc.end != self.lexer.src.len),
            },
        }),

        .invalid => |loc| return diag.ret(error.InvalidBytes, .{
            .line = self.line,
            .loc = loc,
        }),
    };

    const max_triplet_len = comptime @max(
        std.fmt.count("{0d}/{0d}/{0d}", .{std.math.maxInt(IntOneBasedMaybeNeg)}),
        std.fmt.count("{0d}/{0d}/{0d}", .{std.math.minInt(IntOneBasedMaybeNeg)}),
    );
    const TripletBstr = std.BoundedArray(u8, max_triplet_len);

    const triplet_bstr: TripletBstr, //
    const full_loc: Lexer.Token.Loc //
    = src: {
        var triplet_bstr: TripletBstr = .{};

        var str_iter = self.stringIter(first_loc) orelse {
            triplet_bstr.appendSlice(first_loc.getStr(self.lexer.src)) catch
                return diag.ret(error.InvalidReferenceIndexTriplet, .{
                    .line = line,
                    .loc = first_loc,
                });
            break :src .{
                triplet_bstr,
                first_loc,
            };
        };
        var str_peeker = str_iter.toPeeker();

        var state: enum {
            start,
            start_neg,
            skipping_leading_zeroes,
            reading_digits,
        } = .start;
        mainloop: while (true) {
            const segment, _ = str_peeker.next() orelse break :mainloop;
            std.debug.assert(segment.len != 0);

            for (segment) |char| sw: switch (state) {
                .start => switch (char) {
                    '-' => {
                        state = .start_neg;
                        triplet_bstr.append('-') catch {
                            str_iter.consume();
                            return diag.ret(error.InvalidReferenceIndexTriplet, .{
                                .line = line,
                                .loc = str_iter.fullSrcLoc(),
                            });
                        };
                    },
                    else => continue :sw .start_neg, // fallthrough to share code without changing state
                },
                .start_neg => switch (char) {
                    '0' => state = .skipping_leading_zeroes,
                    '/' => {
                        state = .start;
                        triplet_bstr.append('/') catch {
                            str_iter.consume();
                            return diag.ret(error.InvalidReferenceIndexTriplet, .{
                                .line = line,
                                .loc = str_iter.fullSrcLoc(),
                            });
                        };
                    },
                    else => {
                        state = .reading_digits;
                        continue :sw state;
                    },
                },
                .skipping_leading_zeroes => switch (char) {
                    '0' => {},
                    '/' => {
                        state = .start;
                        triplet_bstr.append('0') catch {
                            str_iter.consume();
                            return diag.ret(error.InvalidReferenceIndexTriplet, .{
                                .line = line,
                                .loc = str_iter.fullSrcLoc(),
                            });
                        };
                    },
                    else => {
                        state = .reading_digits;
                        continue :sw state;
                    },
                },
                .reading_digits => {
                    triplet_bstr.append(char) catch {
                        str_iter.consume();
                        return diag.ret(error.InvalidReferenceIndexTriplet, .{
                            .line = line,
                            .loc = str_iter.fullSrcLoc(),
                        });
                    };
                    if (char == '/') {
                        state = .start;
                    }
                },
            };
        }
        str_iter.consume();

        break :src .{ triplet_bstr, str_iter.fullSrcLoc() };
    };
    const triplet_str = triplet_bstr.constSlice();

    var splitter = std.mem.splitScalar(u8, triplet_str, '/');
    const v_str = splitter.first();
    const vt_str_opt = splitter.next();
    const vn_str_opt = splitter.next();
    if (splitter.next() != null or v_str.len == 0) {
        return diag.ret(error.InvalidReferenceIndexTriplet, .{
            .line = line,
            .loc = full_loc,
        });
    }

    const v_raw = std.fmt.parseInt(IntOneBasedMaybeNeg, v_str, 10) catch {
        return diag.ret(error.InvalidReferenceIndexTriplet, .{
            .line = line,
            .loc = full_loc,
        });
    };
    const v = parseRefIdxFromRawInt(v_raw, indices.endOffsetOf(.positions)) orelse {
        return diag.ret(error.InvalidReferenceIndexTriplet, .{
            .line = line,
            .loc = full_loc,
        });
    };

    const vt_raw_opt: ?IntOneBasedMaybeNeg = raw: {
        const vt_str = vt_str_opt orelse break :raw null;
        if (vt_str.len == 0) break :raw null;
        break :raw std.fmt.parseInt(
            IntOneBasedMaybeNeg,
            vt_str,
            10,
        ) catch return diag.ret(error.InvalidReferenceIndexTriplet, .{
            .line = line,
            .loc = full_loc,
        });
    };
    const vt: Ast.RefIndex = vt: {
        const vt_raw = vt_raw_opt orelse break :vt .null;
        break :vt parseRefIdxFromRawInt(
            vt_raw,
            indices.endOffsetOf(.textures),
        ) orelse return diag.ret(error.InvalidReferenceIndexTriplet, .{
            .line = line,
            .loc = full_loc,
        });
    };

    const vn_raw_opt: ?IntOneBasedMaybeNeg = raw: {
        const vn_str = vn_str_opt orelse break :raw null;
        if (vn_str.len == 0) break :raw null;
        break :raw std.fmt.parseInt(
            IntOneBasedMaybeNeg,
            vn_str,
            10,
        ) catch return diag.ret(error.InvalidReferenceIndexTriplet, .{
            .line = line,
            .loc = full_loc,
        });
    };
    const vn: Ast.RefIndex = vn: {
        const vn_raw = vn_raw_opt orelse break :vn .null;
        break :vn parseRefIdxFromRawInt(
            vn_raw,
            indices.endOffsetOf(.textures),
        ) orelse return diag.ret(error.InvalidReferenceIndexTriplet, .{
            .line = line,
            .loc = full_loc,
        });
    };

    const triplet: Ast.FaceTriplet = .{
        .v = v,
        .vt = vt,
        .vn = vn,
    };

    return .{
        .loc = full_loc,
        .value = triplet,
    };
}

const IntOneBasedMaybeNeg = std.math.IntFittingRange(
    -std.math.maxInt(Ast.RefIndex.Int),
    std.math.maxInt(Ast.RefIndex.Int),
);

fn parseRefIdxFromRawInt(
    int_raw: IntOneBasedMaybeNeg,
    end_index: Ast.Node.Index,
) ?Ast.RefIndex {
    const ref_idx: Ast.RefIndex = ref_idx: {
        if (int_raw == 0 or int_raw > std.math.maxInt(Ast.RefIndex.Int) + 1) {
            return null;
        }

        if (int_raw < 0) {
            if (@abs(int_raw) > end_index) return null;
            const end_index_signed: IntOneBasedMaybeNeg = end_index;
            break :ref_idx .from(@intCast(end_index_signed - int_raw));
        } else {
            break :ref_idx .from(@intCast(int_raw - 1));
        }
    };
    if (ref_idx == .null) return null;
    return ref_idx;
}

// -- string iteration & look-ahead stuff -- //

/// If this returns true, that means that `.{ self.nextTok(), self.nextTok() }` would result
/// in `.{ .backslash_nl, .str }`. Returns false otherwise.
/// Does not modify any state.
fn thereIsStringContinuationAfterEscapedNewline(self: *const Parser) bool {
    self.debug.assertPrevTokenKindWas(.str);
    var peeker: Parser = self.*;
    if (peeker.nextTok() != .backslash_nl) return false;
    if (peeker.nextSkipEscapedNewlineChain() != .str) return false;
    return true;
}

/// If necessary, returns an iterator over the disjointed string.
/// Otherwise returns null, indicating to the caller that `first_loc` represents
/// the full string.
fn stringIter(
    self: *Parser,
    first_loc: Lexer.Token.Loc,
) ?StringIter {
    self.debug.assertPrevTokenEql(.{ .str = first_loc });
    if (!self.thereIsStringContinuationAfterEscapedNewline()) return null;
    return .{
        .parser = .{ .ref = self },
        .state = .first_pending,
        .src_loc_joined = first_loc,
    };
}

const StringIter = struct {
    parser: ParserRefOrPeeker,
    /// Did we return the first one yet.
    state: enum {
        first_pending,
        iterating,
        done,
    },
    src_loc_joined: Lexer.Token.Loc,

    const ParserRefOrPeeker = union(enum) {
        ref: *Parser,
        peeker: Parser,
    };

    fn peekLength(self: *const StringIter) usize {
        std.debug.assert(self.state == .first_pending);
        var length: usize = 0;
        var peeker = self.toPeeker();
        while (true) {
            const segment, _ = peeker.next() orelse break;
            length += segment.len;
        }
        return length;
    }

    fn peekEql(self: StringIter, str: []const u8) bool {
        std.debug.assert(self.state == .first_pending);
        var index: usize = 0;
        var peeker = self.toPeeker();
        while (true) {
            const segment, _ = peeker.next() orelse break;
            const str_rest = str[index..];
            if (segment.len > str_rest.len) return false;
            if (!std.mem.eql(u8, segment, str_rest[0..segment.len])) return false;
            index += segment.len;
        }
        std.debug.assert(index <= str.len);
        return index == str.len;
    }

    fn peekToEnum(self: StringIter, comptime Keyword: type) ?Keyword {
        const max_kw_len = comptime blk: {
            var max_kw_len: usize = 0;
            for (@typeInfo(Keyword).@"enum".fields) |field| {
                max_kw_len = @max(max_kw_len, field.name.len);
            }
            break :blk max_kw_len;
        };
        const bstr = self.peekToBoundedArray(max_kw_len) orelse return null;
        return std.meta.stringToEnum(Keyword, bstr.slice());
    }

    fn peekToBoundedArray(
        self: StringIter,
        comptime max_bytes: usize,
    ) ?std.BoundedArray(u8, max_bytes) {
        var peeker = self.toPeeker();
        var bstr: std.BoundedArray(u8, max_bytes) = .{};
        while (true) {
            const segment, _ = peeker.next() orelse break;
            bstr.appendSlice(segment) catch return null;
        }
        return bstr;
    }

    fn toPeeker(self: StringIter) StringIter {
        std.debug.assert(self.state == .first_pending);
        return .{
            .parser = .{ .peeker = switch (self.parser) {
                .ref => |ref| ref.*,
                .peeker => |copy| copy,
            } },
            .state = self.state,
            .src_loc_joined = self.src_loc_joined,
        };
    }

    /// Can only be called after `iter.next()` has returned null.
    fn fullSrcLoc(iter: *const StringIter) Lexer.Token.Loc {
        std.debug.assert(iter.state == .done);
        return iter.src_loc_joined;
    }

    fn consume(iter: *StringIter) void {
        std.debug.assert(iter.state == .first_pending);
        while (iter.next()) |_| {}
    }

    fn next(iter: *StringIter) ?struct { []const u8, Lexer.Token.Loc } {
        const parser = switch (iter.parser) {
            .ref => |ref| ref,
            .peeker => |*peeker| peeker,
        };
        switch (iter.state) {
            .first_pending => {
                return iter.nextFirst();
            },
            .iterating => {
                if (!parser.thereIsStringContinuationAfterEscapedNewline()) {
                    iter.state = .done;
                    return null;
                }

                std.debug.assert(parser.nextTok() == .backslash_nl);
                const next_loc = parser.nextSkipEscapedNewlineChain().str;

                std.debug.assert(iter.src_loc_joined.end < next_loc.start);
                iter.src_loc_joined.end = next_loc.end;
                return .{
                    next_loc.getStr(parser.lexer.src),
                    next_loc,
                };
            },
            .done => {
                unreachable; // illegal to call after receiving null
            },
        }
    }

    fn nextFirst(iter: *StringIter) struct { []const u8, Lexer.Token.Loc } {
        const parser = switch (iter.parser) {
            .ref => |ref| ref,
            .peeker => |*peeker| peeker,
        };

        std.debug.assert(iter.state == .first_pending);
        iter.state = .iterating;

        parser.debug.assertPrevTokenEql(.{ .str = iter.src_loc_joined });
        return .{
            iter.src_loc_joined.getStr(parser.lexer.src),
            iter.src_loc_joined,
        };
    }
};

fn zipUpFaceTriplet(
    gpa: std.mem.Allocator,
    extra: *std.ArrayListUnmanaged(u32),
    triplet: Ast.FaceTriplet,
    want_extra_capacity: enum {
        no_reserve_capacity_for_null_sentinel,
        do_reserve_capacity_for_null_sentinel,
    },
) std.mem.Allocator.Error!void {
    const curr_layout = triplet.layout();
    const extra_capacity: usize = switch (want_extra_capacity) {
        .no_reserve_capacity_for_null_sentinel => 0,
        .do_reserve_capacity_for_null_sentinel => 1,
    };
    try extra.ensureUnusedCapacity(gpa, curr_layout.zippedLength() + extra_capacity);
    switch (curr_layout) {
        .vertex_only => extra.appendAssumeCapacity(
            @bitCast(triplet.v.nonNull().value().?),
        ),
        .only_vt => extra.appendSliceAssumeCapacity(&.{
            @bitCast(triplet.v.nonNull().value().?),
            @bitCast(triplet.vt.nonNull().value().?),
        }),
        .only_vn => extra.appendSliceAssumeCapacity(&.{
            @bitCast(triplet.v.nonNull().value().?),
            @bitCast(triplet.vn.nonNull().value().?),
        }),
        .both_vt_vn => extra.appendSliceAssumeCapacity(&.{
            @bitCast(triplet.v.nonNull().value().?),
            @bitCast(triplet.vt.nonNull().value().?),
            @bitCast(triplet.vn.nonNull().value().?),
        }),
    }
}

fn tryStrToEnum(
    src: []const u8,
    first_loc: Lexer.Token.Loc,
    maybe_str_iter: ?Parser.StringIter,
    comptime E: type,
) ?E {
    var iter = maybe_str_iter orelse {
        const str = first_loc.getStr(src);
        return std.meta.stringToEnum(E, str);
    };
    std.debug.assert(first_loc.eql(iter.src_loc_joined));
    return iter.peekToEnum(E);
}

/// Based off of `std.mem.indexOfScalarPos`, with some tweaks, and with the
/// logic inversed, to seek the first `index` where `slice[index] != value`.
fn indexOfNotScalarPos(comptime T: type, slice: []const T, start_index: usize, value: T) ?usize {
    if (start_index >= slice.len) return null;

    const backend_supports_vectors = switch (@import("builtin").zig_backend) {
        .stage2_llvm, .stage2_c => true,
        else => false,
    };

    var i: usize = start_index;
    if (backend_supports_vectors and
        !std.debug.inValgrind() and // https://github.com/ziglang/zig/issues/17717
        !@inComptime() and
        (@typeInfo(T) == .int or @typeInfo(T) == .float) and std.math.isPowerOfTwo(@bitSizeOf(T)))
    {
        if (std.simd.suggestVectorLength(T)) |block_len| {
            // For Intel Nehalem (2009) and AMD Bulldozer (2012) or later, unaligned loads on aligned data result
            // in the same execution as aligned loads. We ignore older arch's here and don't bother pre-aligning.
            //
            // Use `std.simd.suggestVectorLength(T)` to get the same alignment as used in this function
            // however this usually isn't necessary unless your arch has a performance penalty due to this.
            //
            // This may differ for other arch's. Arm for example costs a cycle when loading across a cache
            // line so explicit alignment prologues may be worth exploration.

            // Unrolling here is ~10% improvement. We can then do one bounds check every 2 blocks
            // instead of one which adds up.
            const Block = @Vector(block_len, T);
            if (i + 2 * block_len < slice.len) {
                const mask: Block = @splat(value);
                while (true) {
                    inline for (0..2) |_| {
                        const block: Block = slice[i..][0..block_len].*;
                        if (std.simd.firstTrue(block != mask)) |first_true_idx| {
                            return i + first_true_idx;
                        }
                        i += block_len;
                    }
                    if (i + 2 * block_len >= slice.len) break;
                }
            }

            // {block_len, block_len / 2} check
            inline for (0..2) |j| {
                const block_x_len = block_len / (1 << j);
                comptime if (block_x_len < 4) break;

                const BlockX = @Vector(block_x_len, T);
                if (i + block_x_len < slice.len) {
                    const mask: BlockX = @splat(value);
                    const block: BlockX = slice[i..][0..block_x_len].*;
                    if (std.simd.firstTrue(block != mask)) |first_true_idx| {
                        return i + first_true_idx;
                    }
                    i += block_x_len;
                }
            }
        }
    }

    return for (slice[i..], i..) |c, j| {
        if (c != value) break j;
    } else null;
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
