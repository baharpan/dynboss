from subprocess import call
import time
import os

#fileIndex = [14,13,12,11,10,9,8,7,6,5,4,3,2,1,0];
fileIndex = [10];
dataDir = "/home/alan/data/bio/ecoli/split/";
Dsk = "/home/alan/exp/bio/dsk/dsk-1.6906/dsk"
#CosmoPack = "/home/alan/exp/bio/cosmoBowe/cosmo-pack"
#CosmoBuild = "/home/alan/exp/bio/cosmoBowe/cosmo-build"
#CosmoBuildDyn = "/home/alan/exp/bio/cosmoBowe/cosmo-build-dyn"
CosmoPack = "../cosmo-pack"
CosmoBuild = "../cosmo-build"
CosmoBuildDyn= "../cosmo-build-dyn"
CosmoAddRemove= "../cosmo-add-delete"
#CosmoQuery = "/home/alan/exp/bio/cosmoBowe/cosmo-query"
OutFile = "./cmpResultsCosmoBowe.txt"

call([Dsk, "temp.fasta", "64"]);
    
call([ CosmoPack, "temp.solid_kmers_binary" ]);
start = time.time();
call([ CosmoBuild, "temp.solid_kmers_binary.packed" ]);
end = time.time();
timeCosmo = end - start;
print timeCosmo, " s";
statinfo = os.stat( "temp.solid_kmers_binary.packed.dbg");
print statinfo.st_size / 1024.0 / 1024.0
start = time.time();
call([ CosmoAddRemove, "temp.solid_kmers_binary.packed" ]);
end = time.time();
timeCosmoDyn = end - start;
print timeCosmoDyn, " s";
statinfo = os.stat( "temp.solid_kmers_binary.packed.dbg");
print statinfo.st_size / 1024.0 / 1024.0

