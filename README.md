# ffpe-filter-legacy

Legacy filter for removing FFPE artifact (C>T/G>A in the XCG context), which was
originally created by Chip Stewart and Lee Lichtenstein.

It is assembled here to serve as a benchmarking baseline.

## Use

### Option 1

Install dependencies (Java 7, Python 2.7, and Matlab Common Runtime 2013a).
Edit `ffpe-filter.sh` to specify path to dependencies.
Then run `ffpe-filter.sh` with required parameters.

### Option 2

Build docker image by

```
docker build -t ffpe-filter-legacy .
```

Then the docker image may be run by

```
docker run --rm -v $(pwd):/root ffpe-filter-legacy /bin/bash
```

Finally, run `ffpe-filter.sh` with required parameters inside the docker shell.


## Remarks

The `oxog-filter-legacy` workflow must be run first, because the several of the
fields annotated by that workflow are required,
e.g. `i_t_ALT_F1R2`, `i_t_ALT_F2R1`, `i_t_REF_F1R2`, `i_t_REF_F2R1`.
