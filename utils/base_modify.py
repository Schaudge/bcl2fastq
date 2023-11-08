#!/usr/bin/env python3
# Written By Schaudge King
# 2022-07-25
# change base content for bcl.bgzf file
# 1. decompress -> 2. modify base -> 3. compress
from contextlib import closing
from Bio import bgzf
position_cluster_dict = {}


def convert_position_to_cluster(local_position_filepath: str):
    cluster_idx = 12  # the prefix header size, should be skipped
    with open(local_position_filepath) as convert_input:
        for position in convert_input:
            position_cluster_dict[position] = cluster_idx
            cluster_idx += 1


def bcl_base_change_workflow(bcl_bgzf_path: str, output_bcl_path: str, block_size: int = 65536):
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
    :return:
    """
    bcl_cluster_base, new_bci_idx = bytearray(), bytearray()

    # 1. first, we need read bcl.bci file into , tile_cluster_idx
    tile_offset_idx, tile_start_offset = [], []
    with open(bcl_bgzf_path + ".bci", "rb") as bci_input:
        all_bci_bytes = bci_input.read()
    assert len(all_bci_bytes) == 200   # each bci has 200 bytes, header size: 8 bytes
    tile_idx = 8
    while tile_idx < 200:
        tile_offset_idx.append(int.from_bytes(all_bci_bytes[tile_idx:tile_idx+2], byteorder="little", signed=False))
        tile_start_offset.append(int.from_bytes(all_bci_bytes[tile_idx+2:tile_idx+8], byteorder="little", signed=False))
        tile_idx += 8

    # 2. read and decompress the raw bcl bgzf bytes into bcl_cluster_base
    tile_offset_matched, tile_cluster_idx = [False] * len(tile_offset_idx), [0] * len(tile_offset_idx)
    bcl_bgzf_offset, bcl_cluster_offset = 0, 0
    with closing(bgzf.BgzfReader(bcl_bgzf_path, "rb")) as bgzf_input:
        bgzf_block_size, bcl_base_block = bgzf_input._block_raw_length, bgzf_input._buffer   # first block was loaded!
        tile_offset_matched[0], tile_cluster_idx[0] = True, 0
        while bgzf_block_size:
            bcl_cluster_base += bcl_base_block
            bcl_bgzf_offset += bgzf_block_size                # record the raw bgzf offset (total read)
            bcl_cluster_offset += len(bcl_base_block)         # record the decompressed bcl offset (total read)
            if bcl_bgzf_offset in tile_start_offset:
                idx = tile_start_offset.index(bcl_bgzf_offset)
                tile_offset_matched[idx], tile_cluster_idx[idx] = True, bcl_cluster_offset
            try:
                bgzf_block_size, bcl_base_block = bgzf._load_bgzf_block(bgzf_input._handle)
            except StopIteration:  # EOF
                bgzf_block_size = 0

    assert all(tile_offset_matched)   # all tile start matched bgzf header

    # 3. change the base in special cluster idx
    # 0x80 -> A (32)  0x81 -> C    0x82 -> G    0x83 -> T
    # bcl_cluster_base[6] = 0x80

    # 4. compress the modified bytearray into bgzf
    compressed_bgzf_size, new_tile_offset = 0, [0]
    written_bcl_size, bcl_cluster_size = block_size, len(bcl_cluster_base)
    with closing(bgzf.BgzfWriter(output_bcl_path, "wb", compresslevel=7)) as bgzf_output:
        while written_bcl_size < bcl_cluster_size:
            compressed_bgzf_size += bgzf_output._write_block(bcl_cluster_base[written_bcl_size-block_size:written_bcl_size])
            written_bcl_size += block_size
            if written_bcl_size in tile_cluster_idx:
                new_tile_offset.append(compressed_bgzf_size)
        # compress the last bcl block
        last_block_start = written_bcl_size - block_size
        compressed_bgzf_size += bgzf_output._write_block(bcl_cluster_base[last_block_start:bcl_cluster_size])

    assert len(new_tile_offset) == len(tile_start_offset)

    # 5. write the new index file bcl.bgzf.bci for bcl.bgzf
    with open(output_bcl_path + ".bci", "wb") as bci_output:
        bci_output.write(all_bci_bytes[:8])
        for tile_offset, block_offset in zip(tile_offset_idx, new_tile_offset):
            virtual_offset = bgzf.make_virtual_offset(block_offset, tile_offset)
            bci_output.write(virtual_offset.to_bytes(8, byteorder="little", signed=False))


if __name__ == "__main__":
    from sys import argv
    bcl_base_change_workflow(argv[1], argv[2])
