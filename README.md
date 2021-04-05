# R2 Estimator

```
Usage: r2-estimator [opts ...] <in.{sav,bcf,vcf.gz}> 

 -h, --help              Print usage
 -o, --output            Path to output file (default: /dev/stdout)
 -O, --output-format     Output file format (sav, bcf, or vcf; default: sav)
 -t, --filter-threshold  Minimum r-squared threshold for output
```

## Build and Install
```shell
# cget available in pip3
cget ignore xz
cget install -f requirements.txt
mkdir build; cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake -DCMAKE_INSTALL_PREFIX=<instal_prefix> ..
make
make install #optional
```
