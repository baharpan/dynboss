from subprocess import call
import time
import os
import sys

Dsk = "/home/alan/exp/bio/dsk/dsk-1.6906/dsk"

CosmoPack = "../cosmo-pack"
CosmoValidate= "../cosmo-validate"

k=27

fnamebase=sys.argv[1];
fileFasta= fnamebase + ".fasta"
call([Dsk, fileFasta , str(k)]);

fileFastaOne= fnamebase + "_1.fasta"
call([Dsk, fileFastaOne , str(k)]);

fileKmer= fnamebase + ".solid_kmers_binary";
fileKmerOne= fnamebase + "_1.solid_kmers_binary";

call([ CosmoPack, fileKmer ]);
call([ CosmoPack, fileKmerOne ]);

filePacked= fnamebase + ".solid_kmers_binary.packed";
filePackedOne = fnamebase + "_1.solid_kmers_binary.packed";

start = time.time();
call([ CosmoValidate, fnamebase ]);
end = time.time();
timeCosmo = end - start;
print timeCosmo, " s";

