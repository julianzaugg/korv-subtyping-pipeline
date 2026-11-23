#!/usr/bin/env python3
"""Extract trim adapters from a CLC Genomics Workbench adapter list (.clc).

A .clc file is a small header, then an LZ4 frame wrapping a Java-serialised
TrimAdapterList. We decode the LZ4 frame (pure Python, no dependencies) and read
the adapter sequences and strands out of the serialised strings.

Usage: extract_clc_adapters.py "ENV adapters.clc" > ENV_primers.fasta
"""
import re
import struct
import sys

LZ4_MAGIC = b"\x04\x22\x4d\x18"


def lz4_block(buf):
    out = bytearray()
    i = 0
    while i < len(buf):
        token = buf[i]
        i += 1
        lit = token >> 4
        if lit == 15:
            while True:
                x = buf[i]
                i += 1
                lit += x
                if x != 255:
                    break
        out += buf[i:i + lit]
        i += lit
        if i >= len(buf):
            break
        offset = buf[i] | (buf[i + 1] << 8)
        i += 2
        mlen = token & 15
        if mlen == 15:
            while True:
                x = buf[i]
                i += 1
                mlen += x
                if x != 255:
                    break
        for _ in range(mlen + 4):
            out.append(out[-offset])
    return bytes(out)


def lz4_frame(data):
    p = data.index(LZ4_MAGIC) + 4
    flg = data[p]
    p += 2  # FLG, BD
    if flg & 0x08:  # content size present
        p += 8
    p += 1  # header checksum
    out = bytearray()
    while True:
        (size,) = struct.unpack("<I", data[p:p + 4])
        p += 4
        if size == 0:
            break
        raw = size & 0x80000000
        size &= 0x7FFFFFFF
        block = data[p:p + size]
        p += size
        if flg & 0x10:  # block checksum
            p += 4
        out += block if raw else lz4_block(block)
    return bytes(out)


def parse_adapters(payload):
    """Return [(name, sequence, strand)] in file order.

    Within the serialised list the entries appear as
    'forward adapters' ... 'Sequence' <seq> ... 'Strand' <Plus|Minus> ... for the
    first adapter, then '<seq>' 'Minus' etc. for the next. The edit history repeats
    some strings, so we take only DNA strings that directly precede a strand word.
    """
    strings = [m.group().decode() for m in re.finditer(rb"[ -~]{3,}", payload)]
    adapters = []
    for idx, s in enumerate(strings):
        if s.startswith(("Plus", "Minus")):
            strand = "Plus" if s.startswith("Plus") else "Minus"
            for j in range(idx - 1, max(idx - 4, -1), -1):
                seq = re.sub(r"^[^ACGTN]+", "", strings[j])
                seq = re.match(r"[ACGTN]{10,}", seq)
                if seq:
                    adapters.append((seq.group(), strand))
                    break
    return adapters


def main(path):
    payload = lz4_frame(open(path, "rb").read())
    adapters = parse_adapters(payload)
    if len(adapters) != 2:
        sys.exit(f"expected 2 adapters, parsed {adapters}")
    names = ["forward_adapter", "reverse_adapter"]
    for name, (seq, strand) in zip(names, adapters):
        print(f">{name} clc_strand={strand}")
        print(seq)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
