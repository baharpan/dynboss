This repository will be updated soon.
# Dynamic BOSS
Dynamic succinct de Bruijn graph

## Dependencies
- SDSL -- [Succinct Data Structure Library](https://github.com/simongog/sdsl-lite)

  To avoid the necessity of editing the Makefile, SDSL may be installed to `/usr/local/`.
- [boost](https://github.com/boostorg/boost)
- [tclap](http://tclap.sourceforge.net/)
- [DSK](https://github.com/GATB/dsk), version 1.6906. The source code for this version is provided in this repository, in `dsk-1.6906`.

Boost and tclap may be installed on a Debian-based system with the following commands.
```
sudo apt-get install libboost-dev
sudo apt-get install libtclap-dev	
```

## How to compile
### Dynamic BOSS
Update the paths on Makefile (unnecessary if dependencies installed on a Debian-based system as described above).
```
make revcomps=0
```
### DSK
```
cd dsk-1.6906

#for Max kmer size 32:
make dsk canon=0

#for larger kmer size like 64:
make dsk k=64 canon=0

```
DSK supports larger `k` than `64`, but currently this is the largest value supported by our implementation.
## Building the data structure
Start with input FASTA file, for example `test/yeast_1.fasta`, which is included. No `N` characters are allowed in the FASTA file.  Run DSK (example with `k=31`):
```
cd test
../dsk-/1.6906/dsk yeast_1.fasta 31
```
This produces output file `yeast_1.solid_kmers_binary`. Next, run `cosmo-pack`:
```
../cosmo-pack yeast_1.solid_kmers_binary
```
which produces output `yeast_1.solid_kmers_binary.packed`
Finally, run `cosmo-build-dyn` to build dynamic BOSS:
```
../cosmo-build-dyn yeast_1.solid_kmers_binary.packed
```
to produce output file `yeast_1.solid_kmers_binary.packed.dbg`. 
## Validation of addition
This test confirms that `add-edge` produces the same graph as the static construction method.
This test constructs BOSS using a static construction method, as well as online from an empty graph
using the dynamic `add_edge` function. The two graphs are then checked to ensure they contain the
same edges and nodes. The script
`test/runAddVerify.bash` runs the test.  Input is a single FASTA file, for example `test/yeast_1.fasta`.
```
cd test
bash runAddVerify.bash yeast_1
```
## Validation of deletion
These tests confirm that deletion is the inverse of addition.
### Test 1
This test loads the constructsd dynamicBOSS, then adds `k`-mers from a FASTA file into the graph with `add-edge` then deletes them with `delete-edge` and confirms that the starting and final graphs are the same.
`test/runDelVerify.bash` runs the test.  Input are the `.dbg` and `.fasta` files. For example `test/yeast_1.solid_kmers_binary.packed.dbg` and `test/yeast_2.fasta`. 
```
cd test
bash runDelVerify.bash yeast_1.solid_kmers_binary.packed.dbg yeast_2.fasta
```

### Test 2
This test deletes random `k`-mers from a starting graph, then adds them back in and verifies that the initial and final graphs are the same.
This test loads the constructsd dynamicBOSS,
then deletes the user-defined number of random `k`-mers with `delete-edge`, adds them back to the graph and confirms that the starting and final graphs are the same.
`test/runDelAdd.bash` runs the test.  Inputs are `.dbg` file and an integer, for example `test/yeast_1.solid_kmers_binary.packed.dbg` and
`1000`.
```
cd test
bash runDelAdd.bash yeast_1.solid_kmers_binary.packed.dbg 1000
```
