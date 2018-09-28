from subprocess import call
import time
import os

#fileIndex = [14,13,12,11,10,9,8,7,6,5,4,3,2,1,0];
fileIndex = [10];
dataDir = "/home/alan/data/bio/ecoli/split/";
Dsk = "/home/alan/exp/bio/dsk/dsk-1.6906/dsk"
CosmoPack = "../cosmo-pack"
CosmoBuild = "../cosmo-build"
CosmoVerify = "../cosmo-verify"
CosmoBuildDyn = "../cosmo-build-dyn"
#CosmoQuery = "/home/alan/exp/bio/cosmoBowe/cosmo-query"
OutFile = "./cmpResultsCosmoBowe.txt"

#call([Dsk, "temp.fasta", "28"]);
call([Dsk, "temp.fasta", "4"]);
    
call([ CosmoPack, "temp.solid_kmers_binary" ]);
call([ CosmoVerify, "temp.solid_kmers_binary.packed" ]);

#call([ CosmoQuery, "temp.solid_kmers_binary.packed.dbg", dataDir + "sra_data" + str( i ) + ".fasta", "-o " + OutFile ]);

# ofObject = open( OutFile, "a" );
# call(["cp", dataDir + "sra_data.fasta", "./temp.fasta" ]);
# call([Dsk, "temp.fasta", "28"]);
    
# start = time.time();
# call([ CosmoPack, "temp.solid_kmers_binary" ]);
# call([ CosmoBuild, "temp.solid_kmers_binary.packed" ]);
# end = time.time();
# timeCosmo = end - start;
# print timeCosmo, " s";
# statinfo = os.stat( "temp.solid_kmers_binary.packed.dbg");
# print statinfo.st_size / 1024.0 / 1024.0
# ofObject.write( "sra_data" + " 27  " + str(timeCosmo) + " " + str(statinfo.st_size / 1024.0 / 1024.0) + " " );
# ofObject.close();
# call([ CosmoQuery, "temp.solid_kmers_binary.packed.dbg", dataDir + "sra_data.fasta", "-o " + OutFile ]);

