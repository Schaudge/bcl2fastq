#!/bin/bash
# Written By Schaudge King
# 2022-07-25

cycle_number="0152"
offset_list="16929" # comma sep for multiple position

if [ ! -e "${cycle_number}.bcl.bgzf.bak" ]; then
    mv ${cycle_number}.bcl.bgzf ${cycle_number}.bcl.bgzf.bak
fi

python3 base_modify.py ${cycle_number} ${offset_list} && bgzip ${cycle_number}.bcl && mv ${cycle_number}.bcl.gz ${cycle_number}.bcl.bgzf

