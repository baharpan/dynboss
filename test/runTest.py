from subprocess import call
import time
import os
import sys

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
#CosmoAddRemove= "../cosmo-delete-add"
#CosmoQuery = "/home/alan/exp/bio/cosmoBowe/cosmo-query"
OutFile = "./cmpResultsCosmoBowe.txt"

fnamebase=sys.argv[1];
fileFasta= fnamebase + ".fasta"
call([Dsk, fileFasta , "4"]);
fileKmer= fnamebase + ".solid_kmers_binary";
filePacked= fnamebase + ".solid_kmers_binary.packed";
call([ CosmoPack, fileKmer ]);
start = time.time();
call([ CosmoBuildDyn, filePacked ]);
end = time.time();
timeCosmo = end - start;
print timeCosmo, " s";
# fileDbg= fnamebase + ".solid_kmers_binary.packed.dbg";
# statinfo = os.stat( fileDbg );
# print statinfo.st_size / 1024.0 / 1024.0
# start = time.time();
# call([ CosmoAddRemove, filePacked ] );
# end = time.time();
# timeCosmoDyn = end - start;
# print timeCosmoDyn, " s";
# statinfo = os.stat( fileDbg );
# print statinfo.st_size / 1024.0 / 1024.0

