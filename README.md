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
Update the paths on Makefile (unnecessary if sdsl-lite is installed to /usr/local/).
```
make
```
## How to validate addition
```
cd dsk-1.6906
make
cd ../test
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
