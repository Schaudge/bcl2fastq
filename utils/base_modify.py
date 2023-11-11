#!/usr/bin/env python3
# Written By Schaudge King
# Init at 2022-07-25
# Change any base or quality content for bcl.bgzf file, and reconstruct the bcl.bgzf.bci index file
# Tested for Illumina MiniSeq Sequencer RTA () output
# Workflow:
# 1. read bci -> 2. decompress bcl (and record flag) -> 3. modify (filter passed) base ->
# 4. compress modified bcl (and index flag) -> 5. write the new bci file
from contextlib import closing
from Bio import bgzf
import struct
position_cluster_dict = {}


def convert_position_to_cluster(local_position_filepath: str):
    cluster_idx = 12  # the prefix header size, should be skipped
    with open(local_position_filepath) as convert_input:
        for position in convert_input:
            position_cluster_dict[position] = cluster_idx
            cluster_idx += 1


def bcl_base_change_workflow(bcl_bgzf_path: str, output_bcl_path: str,
                             block_size: int = 65536, bci_token_size: int = 8):
    """

    The following two byte types were used,
    bytes:     an immutable array of bytes, see help(bytes) for more
    bytearray: a mutable bytearray object, see help(bytearray) for more
    ---------------------------------------------------------------------------------
    with open(bcl_bgzf_path, "rb") as bgzf_input:   # builtins open, not gzip or bgzf
        for each_block in bgzf.BgzfBlocks(bgzf_input):
            print("offset %i, length %i; data start %i, data length %i" % each_block)

    :param bcl_bgzf_path:
    :param output_bcl_path:
    :param block_size: 65536 (2**16) fixed
    :param bci_token_size: 8 bytes (64 bits) fixed
    :return:
    """
    bcl_cluster_base, new_bci_idx = bytearray(), bytearray()

    # 1. first, we need read bcl.bci file into tile_uncompress_padding, next_tile_start_offset
    tile_uncompress_padding, next_tile_start_offset = [], []
    with open(bcl_bgzf_path + ".bci", "rb") as bci_input:
        all_bci_bytes = bci_input.read()
    assert len(all_bci_bytes) == 200   # each bci has 200 bytes, header size: 8 bytes
    tile_idx = 8  # header size
    while tile_idx < 200:
        virtual_offset = struct.unpack("<Q", all_bci_bytes[tile_idx: tile_idx + bci_token_size])[0]  # Q for 64 bits unpack
        start_offset, padding = bgzf.split_virtual_offset(virtual_offset)
        tile_uncompress_padding.append(padding)
        # int.from_bytes(all_bci_bytes[tile_idx: tile_idx + 2], byteorder="little", signed=False)
        next_tile_start_offset.append(start_offset)
        # int.from_bytes(all_bci_bytes[tile_idx + 2: tile_idx + 8], byteorder="little", signed=False)
        tile_idx += bci_token_size

    # 2. read and decompress the raw bcl bgzf bytes into bcl_cluster_base
    tile_offset_matched, tile_cluster_idx = [False] * len(next_tile_start_offset), [-1] * len(next_tile_start_offset)
    bcl_bgzf_offset, bcl_cluster_offset = 0, 0
    with closing(bgzf.BgzfReader(bcl_bgzf_path, "rb")) as bgzf_input:
        bgzf_block_size, bcl_base_block = bgzf_input._block_raw_length, bgzf_input._buffer   # first block was loaded!
        tile_offset_matched[0] = True
        while bgzf_block_size:
            if bcl_bgzf_offset in next_tile_start_offset:
                idx = next_tile_start_offset.index(bcl_bgzf_offset)
                tile_offset_matched[idx], tile_cluster_idx[idx] = True, bcl_cluster_offset
            bcl_cluster_base += bcl_base_block
            bcl_bgzf_offset += bgzf_block_size                # record the raw bcl bgzf offset (total read)
            bcl_cluster_offset += len(bcl_base_block)         # record the decompressed bcl offset (total read)
            try:
                bgzf_block_size, bcl_base_block = bgzf._load_bgzf_block(bgzf_input._handle)
            except StopIteration:  # EOF
                bgzf_block_size = 0

    assert all(tile_offset_matched)   # all bgzf.bci tile start matched bcl.bgzf block

    # 3. change the base in special cluster idx
    # 0x80 -> A (32)  0x81 -> C    0x82 -> G    0x83 -> T
    # bcl_cluster_base[6] = 0x80

    # 4. compress the modified bytearray into bgzf
    compressed_bgzf_size, new_tile_offset = 0, [0]
    written_bcl_size, bcl_cluster_size = block_size, len(bcl_cluster_base)
    with closing(bgzf.BgzfWriter(output_bcl_path, "wb", compresslevel=7)) as bgzf_output:
        while written_bcl_size < bcl_cluster_size:
            compressed_bgzf_size += bgzf_output._write_block(bcl_cluster_base[written_bcl_size-block_size:written_bcl_size])
            if written_bcl_size in tile_cluster_idx:
                new_tile_offset.append(compressed_bgzf_size)
            written_bcl_size += block_size   # ... next decompress block end index
        # compress the last bcl block
        last_block_start = written_bcl_size - block_size
        compressed_bgzf_size += bgzf_output._write_block(bcl_cluster_base[last_block_start:bcl_cluster_size])

    assert len(new_tile_offset) == len(tile_uncompress_padding)

    # 5. write the new index file bcl.bgzf.bci for bcl.bgzf
    with open(output_bcl_path + ".bci", "wb") as bci_output:
        bci_output.write(all_bci_bytes[:bci_token_size])
        for tile_offset, block_offset in zip(tile_uncompress_padding, new_tile_offset):
            virtual_offset = bgzf.make_virtual_offset(block_offset, tile_offset)
            # bci_output.write(struct.pack("<Q", virtual_offset))
            bci_output.write(virtual_offset.to_bytes(bci_token_size, byteorder="little", signed=False))


if __name__ == "__main__":
    from sys import argv
    bcl_base_change_workflow(argv[1], argv[2])
