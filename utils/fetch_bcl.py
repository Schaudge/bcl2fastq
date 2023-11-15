#!/usr/bin/env python3
# Written By Schaudge King
# 2023-11-03
from os import path, makedirs, listdir
position_cluster_dict = {}


def convert_position_to_cluster(local_position_filepath: str):
    cluster_idx = 4  # the prefix bcl header size, should be skipped
    with open(local_position_filepath) as convert_input:
        for position in convert_input:
            position_cluster_dict[position.strip()] = cluster_idx
            cluster_idx += 1


def fetch_bcl_position_cluster(reads_name_file: str, input_bcl_dir: str, output_bcl_dir: str):
    """
    fetch bcl bytes from special bcl position (cluster index)
    :param reads_name_file:
    :param input_bcl_dir:
    :param output_bcl_dir:
    :return:
    """
    position_cluster_idx = []
    with open(reads_name_file) as reads_name_input:
        for read_name in reads_name_input:
            if read_name[0].isdigit():
                position_cluster_idx.append(position_cluster_dict[read_name.strip()])
            else:
                __, __, __, lane, tile, pos_x, pos_y = read_name.strip().split(":")
                cluster_key = lane + ":" + tile + ":" + pos_x + ":" + pos_y
                cluster_idx = position_cluster_dict[cluster_key]
                if cluster_idx not in position_cluster_idx:
                    position_cluster_idx.append(cluster_idx)
    position_cluster_idx.sort()

    # fetch each cluster for each cycle bcl data
    for bcl in listdir(input_bcl_dir):
        if bcl.endswith(".bcl"):
            output_top_dir = output_bcl_dir + "/" + path.basename(reads_name_file).split(".")[0] + "/"
            makedirs(output_top_dir, exist_ok=True)
            with open(input_bcl_dir + "/" + bcl, "rb") as input_bcl, open(output_top_dir + bcl, "wb") as bcl_output:
                all_bcl_bytes = input_bcl.read()
                fetched_bcl_bytes = bytearray()
                for idx in position_cluster_idx:
                    fetched_bcl_bytes += bytes(all_bcl_bytes[idx: idx + 1])
                bcl_output.write(fetched_bcl_bytes)


if __name__ == "__main__":
    from sys import argv
    local_convert_table, fetch_reads_name, raw_bcl_dir, new_bcl_dir = argv[1], argv[2], argv[3], argv[4]
    convert_position_to_cluster(local_convert_table)
    fetch_bcl_position_cluster(fetch_reads_name, raw_bcl_dir, new_bcl_dir)

