# modify the base and quality context in bcl file

## 1. python3 and bgzip install

## 2. copy scripts `bcl_correct.sh` and `base_modify.py` into the bcl (e.g. Data/Intensities/BaseCalls/L001) directory

## 3. setting the bcl_correct.sh arguments for execution (see bcl_correct.sh for details)
```
sh bcl_correct.sh cycle_number offset_list
```

## 4. rerun bcl2fastq for checking
