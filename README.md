# Dynamic BOSS
Dynamic succinct de Bruijn graph

## How to compile
```
#change the paths on Makefile
make DYNAMIC=1
```
## How to run
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
