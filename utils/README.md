# modify the base and quality context in bcl file

## 1. bgzip all bcl.bgzf blocks into bcl cluster bytes.

## 2. fetch (construct) the positive sites cycle bcl database (one dir to one site), see fetch_bcl.py for reference.

## 3. run the base modify workflow for one batch from configuration.
```
python3 base_modify.py sites_conf_path
```

## 4. rerun bcl2fastq for checking
