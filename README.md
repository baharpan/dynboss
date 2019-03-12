# Dynamic BOSS
Dynamic succinct de Bruijn graph

## Dependencies
- SDSL -- [Succinct Data Structure Library](https://github.com/simongog/sdsl-lite)

  To avoid the necessity of editing the Makefile, SDSL may be installed to `/usr/local/`.
- [boost](https://github.com/boostorg/boost)
- [tclap](http://tclap.sourceforge.net/)

Boost and tclap may be installed on a Debian-based system with the following commands.
```
sudo apt-get install libboost-dev
sudo apt-get install libtclap-dev	
```

## How to compile
### Dynamic BOSS
Update the paths on Makefile (unnecessary if dependencies installed on a Debian-based system as described above).
```
make
```
### DSK
```
cd dsk-1.6906
make k=64
```
## Building the data structure
Start with input FASTA file, for example `test/yeast_1.fasta`, which is included. No `N` characters are allowed in the FASTA file.  Run DSK (example with `k=40`):
```
cd test
../dsk-/1.6906/dsk yeast_1.fasta 40
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
## How to validate addition
This test constructs BOSS using a static construction method, as well as online from an empty graph
using the dynamic `add_edge` function. The two graphs are then checked to ensure they contain the
same edges and nodes. The script
`test/runAddVerify.bash` runs the test.  Input is a single FASTA file, for example `test/yeast_1.fasta`.
```
cd test
bash runAddVerify.bash yeast_1
```
## How to test
```
#count the k-mers
mkdir -p kmc_temp
ls -1 --color=no *.fasta |xargs -l -i echo "~/kmc -b -fq -k32 -ci0 -cs250 {} {}.kmc kmc_temp" >kmercount.sh
source kmercount.sh
ls -1 --color=no *.fasta |xargs -l -i echo "~/kmc_tools sort {}.kmc {}.kmc.sorted " >kmercountsort.sh
source kmercountsort.sh
ls -1 --color=no *.fasta |xargs -l -i echo "{}.kmc.sorted" > filtered_kmc2_list

#build BOSS
./cosmo-pack -k filtered_kmc2_list
#build Dynamic BOSS
./cosmo-add-delete filtered_kmc2_list.packed
```
