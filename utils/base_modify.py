#!/usr/bin/env python3
# Written By Schaudge King
# 2022-07-25
# change base for bcl file
import gzip


base_convert_dict = {"A": 148, "C": 150, "G": 149, "T": 151, "N": 0}


def bcl_base_change(cycle_num, offset_list, base_list):
    bak_bcl_path = "{}.bcl.bgzf.bak".format(cycle_num)
    std_bcl_path = "{}.bcl".format(cycle_num)
    with gzip.open(bak_bcl_path, "rb") as bcl_input, open(std_bcl_path, "wb") as output:
        bcl_byte_stream = bcl_input.read()
        bcl_char_buffer = bytearray(bcl_byte_stream)
        for offset, base in zip(offset_list, base_list):
            print("offset: {} with byte {}".format(offset, bcl_char_buffer[offset]))
            bcl_char_buffer[offset] = base_convert_dict[base]
        output.write(bcl_char_buffer)


if __name__ == "__main__":
    from sys import argv
    offset_list = []
    base_list = []
    if len(argv) > 3:
        offset_list = [int(offset) for offset in argv[2].split(",")]
        base_list = argv[3].split(",")
        if len(offset_list) != len(base_list):
            print("make sure the offset_list and base_list have same size!")
            exit(-1)
    elif len(argv) > 2:
        offset_list = [int(offset) for offset in argv[2].split(",")]
        base_list = ["N"] * len(offset_list)
    bcl_base_change(argv[1], offset_list, base_list)

