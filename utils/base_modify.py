#!/usr/bin/env python3
# Written By Schaudge King
# Init at 2022-07-25
# Change any base or quality content for bcl.bgzf file, and reconstruct the bcl.bgzf.bci index file
# Tested for Illumina MiniSeq Sequencer RTA () output
# Workflow:
# 1. read bci -> 2. decompress bcl (and record flag) -> 3. modify (filter passed) base ->
# 4. compress modified bcl (and index flag) -> 5. write the new bci file
from os import path
from random import sample
from contextlib import closing
from Bio import bgzf
from json import load
import struct


def read_correct_sites_conf(batch_dir_path: str, sites_correct_conf: str, cycle_bcl_name: str,
                            first_cycle: bool = False, adapter_base: bool = False):
    # 1. read sample <-> sites correction configration for one batch!
    with open(sites_correct_conf) as conf_input:
        sites_conf_dict = load(conf_input)

    # 2. read all Undetermined FILTER PASS cluster index for one batch
    all_idx_pool, sampled_idx_pool = [], []
    with open(batch_dir_path + "/UN") as unknown_idx_input:
        for _cluster_idx in unknown_idx_input:
            all_idx_pool.append(int(_cluster_idx.strip()))

    # 3. build idx <-> base quality pair, and output the random setting for first cycle data ...
    correct_sites_dict, filter_sites_dict = {}, {}
    for sample_id, (site, reads_count, cluster_path) in sites_conf_dict.items():
        cluster_file_path = batch_dir_path + "/" + cluster_path
        fixed_random_setting = batch_dir_path + "/" + sample_id + "_" + site
        select_sites_file = ("adapter/" + sample_id + "/" + cycle_bcl_name) if adapter_base else (site + "/" + cycle_bcl_name)
        if cluster_path == "UN":
            with open(select_sites_file, "rb") as bcl_input:
                all_positive_sites = bytearray(bcl_input.read())
            if first_cycle:
                sampled_cluster_idx = sample(range(0, len(all_idx_pool)), int(reads_count))
                random_base_idx = sample(range(0, len(all_positive_sites)), int(reads_count))
                with open(fixed_random_setting, "w") as random_setting_output:
                    for _idx1, _idx2 in zip(sampled_cluster_idx, random_base_idx):
                        while _idx1 in sampled_idx_pool:
                            _idx1 += 1
                        correct_sites_dict[all_idx_pool[_idx1]] = all_positive_sites[_idx2: _idx2 + 1]
                        sampled_idx_pool.append(_idx1)
                        random_setting_output.write(str(_idx1) + "\t" + str(_idx2) + "\n")
            elif path.exists(fixed_random_setting):
                with open(fixed_random_setting) as idx_base_input:
                    for line in idx_base_input:
                        for _idx1, _idx2 in line.strip().split("\t"):
                            fixed_idx1, fixed_idx2 = int(_idx1), int(_idx2)
                            correct_sites_dict[all_idx_pool[fixed_idx1]] = all_positive_sites[fixed_idx2: fixed_idx2 + 1]
            else:
                raise Exception("break UN run ---")
        elif cluster_path == "FILTER":
            with open(cluster_file_path) as idx_input:
                for line in idx_input:
                    filter_sites_dict[int(line.strip())] = True
        elif path.exists(cluster_file_path):
            inplace_cluster_idx = []
            with open(cluster_file_path) as idx_input:
                for line in idx_input:
                    inplace_cluster_idx.append(int(line.strip()))
            with open(select_sites_file, "rb") as bcl_input:
                all_positive_sites = bytearray(bcl_input.read())
            random_base_idx = sample(range(0, len(all_positive_sites)), int(reads_count))
            if first_cycle:
                with open(fixed_random_setting, "w") as random_setting_output:
                    for _idx1, _idx2 in zip(inplace_cluster_idx, random_base_idx):
                        correct_sites_dict[_idx1] = all_positive_sites[_idx2: _idx2 + 1]
                        random_setting_output.write(str(_idx1) + "\t" + str(_idx2) + "\n")
            elif path.exists(fixed_random_setting):
                with open(fixed_random_setting) as idx_base_input:
                    for line in idx_base_input:
                        for _idx1, _idx2 in line.strip().split("\t"):
                            fixed_idx1, fixed_idx2 = int(_idx1), int(_idx2)
                            correct_sites_dict[fixed_idx1] = all_positive_sites[fixed_idx2: fixed_idx2 + 1]
            else:
                raise Exception("break inplace run ---")
    return correct_sites_dict, filter_sites_dict


def bcl_base_change_workflow(bcl_bgzf_path: str, output_bcl_path: str, sites_conf_path: str,
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
    :param sites_conf_path:
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
    first_sequencing_cycle = True if bcl_bgzf_path.endswith("0001.bcl.bgzf") else False
    adapter_cycle_flag = True if path.basename(bcl_bgzf_path) in ["0152.bcl.bgzf", "0153.bcl.bgzf", "0154.bcl.bgzf",
                                                                  "0155.bcl.bgzf", "0156.bcl.bgzf", "0157.bcl.bgzf",
                                                                  "0158.bcl.bgzf", "0159.bcl.bgzf"] else False
    inplace_sites, filter_sites = read_correct_sites_conf(bcl_bgzf_path, sites_conf_path, path.basename(bcl_bgzf_path),
                                                          first_sequencing_cycle, adapter_cycle_flag)
    for cluster_idx, quality_base_byte in inplace_sites.items():
        bcl_cluster_base[cluster_idx] = quality_base_byte
    if first_sequencing_cycle:  # change the filter flag
        input_filter_file = path.dirname(bcl_bgzf_path) + "/s_1.filter"
        output_filter_file = path.dirname(output_bcl_path) + "/s_1.filter"
        with open(input_filter_file, "rb") as filter_input, open(output_filter_file, "wb") as filter_output:
            filter_byte_stream = filter_input.read()
            filter_char_buffer = bytearray(filter_byte_stream)
        for cluster_idx in filter_sites:
            filter_char_buffer[cluster_idx] = 0x0
        filter_output.write(filter_char_buffer)

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
    input_bcl_file, output_bcl_file, correct_conf_file = argv[1], argv[2], argv[3]
    bcl_base_change_workflow(input_bcl_file, output_bcl_file, correct_conf_file)
