#!/usr/bin/env python3
# Written By Schaudge King
# 2023-11-03
# change filter flag state for s_1.filter file
# @TPNB500450:928:HTWJWBGXT:1:11101:3581:1144 1:N:0:AATGTTGC
#                                               |
#                                    0: filter, 1: keep
position_cluster_dict = {}


def convert_position_to_cluster(local_position_filepath: str):
    cluster_idx = 12  # the prefix header size, should be skipped
    with open(local_position_filepath) as convert_input:
        for position in convert_input:
            position_cluster_dict[position] = cluster_idx
            cluster_idx += 1


def cluster_filter_change(reads_name_file: str, input_filter_file: str, output_filter_file: str):
    """
    change the filter flag to 0 for inputs reads name
    :param reads_name_file:
    :param input_filter_file:
    :param output_filter_file:
    :return:
    """
    cluster_idx_list = []
    with open(reads_name_file) as reads_name_input:
        for read_name in reads_name_input:
            if read_name[0].isdigit():
                cluster_idx_list.append(position_cluster_dict[read_name])
            else:
                __, __, __, lane, tile, pos_x, pos_y = read_name.strip().split(":")
                cluster_key = lane + ":" + tile + ":" + pos_x + ":" + pos_y
                cluster_idx = position_cluster_dict[cluster_key]
                if cluster_idx not in cluster_idx_list:
                    cluster_idx_list.append(cluster_idx)

    with open(input_filter_file, "rb") as filter_input, open(output_filter_file, "wb") as filter_output:
        filter_byte_stream = filter_input.read()
        filter_char_buffer = bytearray(filter_byte_stream)
        for cluster_idx in cluster_idx_list:
            print(cluster_idx_list[cluster_idx - 1], cluster_idx_list[cluster_idx], cluster_idx_list[cluster_idx + 1])
            filter_char_buffer[cluster_idx] = 0x0
        filter_output.write(filter_char_buffer)


if __name__ == "__main__":
    from sys import argv

    local_convert_table, filter_reads_name, raw_filter_filename, new_filter_filename = argv[1], argv[2], argv[3], argv[4]
    convert_position_to_cluster(local_convert_table)
    cluster_filter_change(filter_reads_name, raw_filter_filename, new_filter_filename)

