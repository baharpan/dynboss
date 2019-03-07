# Dynamic BOSS
Dynamic succinct de Bruijn graph

## Dependencies
. sdsl-lite
. DSK
The remaining dependencies may be installed on a Debian-based system with the following commands.
```
	sudo apt-get install libboost-dev
	sudo apt-get install libtclap-dev
	
```

## How to compile
```
#change the paths on Makefile (although if sdsl-lite is installed to /usr/local/, this shouldn't be necessary.
make cosmo-add-verify
```
## How to validate addition
Update `test/runAddVerify.bash` to point to your binary for DSK.
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
