# DynamicBOSS
Succinct Dynamic de Bruijn Graph

## Dependencies
- SDSL -- [Succinct Data Structure Library](https://github.com/simongog/sdsl-lite)
- [boost](https://github.com/boostorg/boost)
- [tclap](http://tclap.sourceforge.net/)
- [DSK](https://github.com/GATB/dsk), version 1.6906. The source code for this version is provided in this repository, in `dsk-1.6906`.
- [DYNAMIC](https://github.com/xxsds/DYNAMIC) and [hopscotch-map](https://github.com/Tessil/hopscotch-map), both added as submodules   in this repository. 

## Installation
```
git clone --recurse-submodules https://github.com/baharpan/dynboss.git 
cd dynboss

#compile DSK
cd dsk-1.6906
#for Max kmer size 32:
make dsk canon=0
#for larger kmer size like 64:
make dsk k=64 canon=0

#compile DynamicBOSS
cd src
#Update the paths on Makefile (sdsl and tclap installed in the folder 3rd_party_inst, and boost in 3rd_party_inst/boost)
make revcomps=0
```
## Input
Input is fasta file. No `N` characters are allowed in the fasta file.  
## kmer counting and packing the edges
```
dsk-1.6906/dsk <.fasta file> <k value>
bin/cosmo-pack <.solid_kmers_binary file>
```
## DynamicBOSS Interface
```
bin/dynamicBOSS   build   -p <.packed file>
bin/dynamicBOSS   add     -g <graph file> -s <kmer file>
bin/dynamicBOSS   delete  -g <graph file> -s <kmer file>
bin/dynamicBOSS   query   -g <graph file> -s <kmer file>
```
`build` builds the graph (`.dbg` file)
`add` counts the distinct kmers in the`<kmer file>`, adds them to the `<graph file>` and writes the updated graph with extesion `.updated`.
`delete` counts the distinct kmers in the`<kmer file>`, deletes them from the `<graph file>` and writes the updated graph with extesion `.updated`.
`query` counts the distinct kmers in the`<kmer file>`, query them in the `<graph file>` and writes the results in file queryResults.tsv.
`<graph file>` must have extensions `.dbg` or `.updated`.
`<kmer file>` must be fasta file with extensions `.fasta` or `.fa`
## Compelete Example
```
cd test
#count kmers with k=31, and output ".solid_kmers_binary" file
../dsk-/1.6906/dsk  yeast_1.fasta 31

#create file "yeast_1.solid_kmers_binary.packed" needed for building DynamicBOSS
../bin/cosmo-pack   yeast_1.solid_kmers_binary

#build the DynamicBOSS called "yeast_1.solid_kmers_binary.packed.dbg"
../bin/dynamicBOSS  build  -p   yeast_1.solid_kmers_binary.packed

#add the distinct kmers in fasta file "yeast_2.fasta" to graph "yeast_1.solid_kmers_binary.packed.dbg" and writes the resulting graph "yeast_1.solid_kmers_binary.packed.dbg.updated"
../bin/dynamicBOSS  add    -g   yeast_1.solid_kmers_binary.packed.dbg -s yeast_2.fasta

#delete the distinct kmers in fasta file "yeast_2.fasta" from graph "yeast_1.solid_kmers_binary.packed.dbg" and writes the resulting graph "yeast_1.solid_kmers_binary.packed.dbg.updated"
../bin/dynamicBOSS  delete -g   yeast_1.solid_kmers_binary.packed.dbg -s yeast_2.fasta

#query the distinct kmers in fasta file "yeast_2.fasta" in graph "yeast_1.solid_kmers_binary.packed.dbg" and write the results in queryResults.tsv
../bin/dynamicBOSS  query  -g   yeast_1.solid_kmers_binary.packed.dbg -s yeast_2.fasta
```

## Validations
To validate the addition and deletion in DynamicBOSS, try below `.bash` files in folder `test`

```
bash runAddVerify.bash <.dbg file> <.fasta file that DynamicBOSS is built on>
bash runDelAdd.bash <.dbg file> <number of kmers to delete and add>
bash runDelVerify.bash <.dbg file> <.fasta file>
```
`runAddVerify.bash` confirms that adding kmers to an empty graph builds the same graph as `dynamicBOSS build`.
`runDelAdd.bash` and `runDelVerify.bash` confirm that deletion is the inverse of addition.

