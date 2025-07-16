const std = @import("std");

const Lexer = @This();
src: []const u8,
index: usize,

pub fn init(src: []const u8) Lexer {
    return .{
        .src = src,
        .index = 0,
    };
}

pub const Token = union(Kind) {
    backslash_invalid: Loc,
    invalid: Loc,
    whitespace: Loc,
    str: Loc,

    backslash_nl: Loc,

    /// Ends with `'\r'` or `'\n'`.
    comment_nl_one: Loc,
    /// Ends with `'\r\n'`.
    comment_nl_two: Loc,
    /// Is one of: `'\r'`, `'\n'`, `\r\n`.
    nl: Loc,

    comment_eof: Loc,
    eof: usize,

    pub fn eql(a: Token, b: Token) bool {
        if (a != b.getKind()) return false;
        return switch (a) {
            inline else => |payload_a, tag| blk: {
                const payload_b = @field(b, @tagName(tag));
                break :blk switch (@TypeOf(payload_a)) {
                    Loc => payload_a.eql(payload_b),
                    usize => payload_a == payload_b,
                    else => comptime unreachable,
                };
            },
        };
    }

    pub fn getKind(self: Token) Kind {
        return self;
    }

    pub const Kind = enum {
        backslash_invalid,
        invalid,
        whitespace,
        str,

        backslash_nl,

        comment_nl_one,
        comment_nl_two,
        nl,

        comment_eof,
        eof,
    };

    pub fn getLoc(self: Token) Loc {
        return switch (self) {
            .invalid,
            .backslash_invalid,
            .whitespace,
            .str,

            .comment_eof,
            .comment_nl_one,
            .comment_nl_two,

            .backslash_nl,
            .nl,
            => |loc| loc,

            .eof => |start| .initAbs(start, comptime 1),
        };
    }

    pub const Loc = extern struct {
        start: usize,
        end: usize,

        pub fn initRel(start: usize, len: usize) Loc {
            return .initAbs(start, start + len);
        }

        pub fn initAbs(start: usize, end: usize) Loc {
            std.debug.assert(start <= end);
            return .{
                .start = start,
                .end = end,
            };
        }

        pub fn getStr(loc: Loc, str: anytype) @TypeOf(str) {
            return str[loc.start..loc.end];
        }

        pub fn length(loc: Loc) usize {
            return loc.end - loc.start;
        }

        pub fn eql(a: Loc, b: Loc) bool {
            return a.start == b.start and a.end == b.end;
        }
    };

    pub fn format(
        self: Token,
        writer: *std.Io.Writer,
    ) std.io.Writer.Error!void {
        try std.zon.stringify.serialize(self, .{ .whitespace = true }, writer);
    }
};

pub fn next(self: *Lexer) Token {
    const src = self.src;

    if (self.index == src.len) {
        return .{ .eof = self.index };
    }

    switch (typedChar(src[self.index])) {
        .space, .tab => {
            const start = self.index;
            self.index += 1;
            self.index = std.mem.indexOfNonePos(u8, src, self.index, " \t") orelse src.len;
            return .{ .whitespace = .initAbs(start, self.index) };
        },

        .lf, .cr => {
            const start = self.index;
            _ = self.processNewline();
            return .{ .nl = .initAbs(start, self.index) };
        },

        .backslash => {
            const start = self.index;
            self.index += 1;

            if (self.index == src.len) {
                return .{ .backslash_invalid = .initAbs(start, self.index) };
            }
            switch (typedChar(src[self.index])) {
                else => {
                    return .{ .backslash_invalid = .initAbs(start, self.index) };
                },
                .lf, .cr => {
                    _ = self.processNewline();
                    return .{ .backslash_nl = .initAbs(start, self.index) };
                },
            }
        },

        .hashtag => {
            const start = self.index;
            self.index += 1;
            while (true) {
                self.index = std.mem.indexOfAnyPos(u8, src, self.index, "\r\n" ++ "\\") orelse src.len;

                if (self.index == src.len) {
                    return .{ .comment_eof = .initAbs(start, self.index) };
                }

                if (src[self.index] == '\\' and self.index + 1 != src.len) {
                    switch (typedChar(src[self.index + 1])) {
                        .lf, .cr => {
                            self.index += 1;
                            _ = self.processNewline();
                            continue;
                        },
                        else => {},
                    }
                }

                return if (self.processNewline())
                    .{ .comment_nl_two = .initAbs(start, self.index) }
                else
                    .{ .comment_nl_one = .initAbs(start, self.index) };
            }
        },

        _ => |char_typed| {
            const start = self.index;
            self.index += 1;

            if (!char_typed.isPrint()) {
                self.index = for (
                    src[self.index..],
                    self.index..,
                ) |next_char, i| {
                    if (std.ascii.isPrint(next_char)) break i;
                } else src.len;
                return .{ .invalid = .initAbs(start, self.index) };
            } else {
                self.index = for (
                    src[self.index..],
                    self.index..,
                ) |next_char, i| {
                    const next_char_typed = typedChar(next_char);
                    switch (next_char_typed) {
                        .hashtag,
                        .space,
                        .tab,
                        .cr,
                        .lf,
                        .backslash,
                        => break i,
                        else => {},
                    }
                    if (!next_char_typed.isPrint()) break i;
                } else src.len;
                return .{ .str = .initAbs(start, self.index) };
            }
        },
    }
}

/// Caller should save `const start = self.index;` before calling this.
/// After calling this, the caller can have `const end = self.index;`,
/// with `src[start..end]` forming the newline, which is either cr, lf,
/// or crlf.
/// Returns true for crlf, false for cr or lf.
fn processNewline(self: *Lexer) bool {
    const src = self.src;
    switch (typedChar(src[self.index])) {
        else => unreachable, // should only call when it's one of lf or cr
        .lf => {
            self.index += 1;
            return false;
        },
        .cr => {
            self.index += 1;
            if (self.index == src.len) return false;
            if (src[self.index] != '\n') return false;
            self.index += 1;
            return true;
        },
    }
}

fn typedChar(byte: u8) TypedChar {
    return @enumFromInt(byte);
}

const TypedChar = enum(u8) {
    hashtag = '#',
    space = ' ',
    tab = '\t',
    cr = '\r',
    lf = '\n',
    backslash = '\\',
    _,

    fn isPrint(char: TypedChar) bool {
        const result = std.ascii.isPrint(@intFromEnum(char));
        switch (char) {
            .hashtag,
            .space,
            .tab,
            .cr,
            .lf,
            .backslash,
            => std.debug.assert(result),
            _ => {},
        }
        return result;
    }
};

test Lexer {
    try testLexer("", &.{});
    try testLexer("\n", &.{.{ .nl, "\n" }}); // LF line ending (linux)
    try testLexer("\r", &.{.{ .nl, "\r" }}); // CR line ending (osx)
    try testLexer("\r\n", &.{.{ .nl, "\r\n" }}); // CRLF line ending (windows)

    try testLexer("\n\r", &.{ .{ .nl, "\n" }, .{ .nl, "\r" } });
    try testLexer("\r\r", &.{ .{ .nl, "\r" }, .{ .nl, "\r" } });
    try testLexer("\n\n", &.{ .{ .nl, "\n" }, .{ .nl, "\n" } });
    try testLexer("\n\r\n", &.{ .{ .nl, "\n" }, .{ .nl, "\r\n" } });
    try testLexer("\r\n\n", &.{ .{ .nl, "\r\n" }, .{ .nl, "\n" } });
    try testLexer("\r\n\r", &.{ .{ .nl, "\r\n" }, .{ .nl, "\r" } });

    try testLexer("\x00", &.{.{ .invalid, "\x00" }});
    try testLexer("\x00\x01", &.{.{ .invalid, "\x00\x01" }});
    try testLexer("\x00 \x01", &.{ .{ .invalid, "\x00" }, .{ .whitespace, " " }, .{ .invalid, "\x01" } });

    try testLexer("fizz # foo bar", &.{
        .{ .str, "fizz" },
        .{ .whitespace, " " },
        .{ .comment_eof, "# foo bar" },
    });
    try testLexer("fizz # foo bar\n", &.{
        .{ .str, "fizz" },
        .{ .whitespace, " " },
        .{ .comment_nl_one, "# foo bar\n" },
    });
    try testLexer("v 1.0 -1.0 1.0", &.{
        .{ .str, "v" },
        .{ .whitespace, " " },
        .{ .str, "1.0" },
        .{ .whitespace, " " },
        .{ .str, "-1.0" },
        .{ .whitespace, " " },
        .{ .str, "1.0" },
    });
    try testLexer(
        \\v 1.0 -1.0 1.0
        \\v 1.0 -1.0 1.0 0.5
    , &.{
        .{ .str, "v" },
        .{ .whitespace, " " },
        .{ .str, "1.0" },
        .{ .whitespace, " " },
        .{ .str, "-1.0" },
        .{ .whitespace, " " },
        .{ .str, "1.0" },
        .{ .nl, "\n" },
        .{ .str, "v" },
        .{ .whitespace, " " },
        .{ .str, "1.0" },
        .{ .whitespace, " " },
        .{ .str, "-1.0" },
        .{ .whitespace, " " },
        .{ .str, "1.0" },
        .{ .whitespace, " " },
        .{ .str, "0.5" },
    });

    try testLexer("foo\\\nbar", &.{
        .{ .str, "foo" },
        .{ .backslash_nl, "\\\n" },
        .{ .str, "bar" },
    });
}

const TestToken = struct { Lexer.Token.Kind, []const u8 };
fn testLexer(
    src: []const u8,
    expected_tokens: []const TestToken,
) !void {
    const gpa = std.testing.allocator;

    var actual_tokens: std.ArrayListUnmanaged(TestToken) = .empty;
    defer actual_tokens.deinit(gpa);

    var lexer: Lexer = .init(src);
    while (true) {
        const tok = lexer.next();
        if (tok == .eof) break;
        try actual_tokens.append(gpa, .{ tok, tok.getLoc().getStr(src) });
    }

    try std.testing.expectEqualDeep(expected_tokens, actual_tokens.items);
}
