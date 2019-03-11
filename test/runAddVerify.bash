
dsk="../dsk-1.6906/dsk"
CosmoPack="../cosmo-pack"
CosmoAddVerify="../cosmo-add-verify"
k=64

#OutFile = "./cmpResultsCosmoBowe.txt"

if [ $# -lt 1 ]; then
    echo "Usage: $0 <fastafilebase>"
    exit
fi

fnamebase="$1"
fileFasta="$1.fasta"

cmd="$dsk $fileFasta $k"
echo $cmd
$cmd

fileKmer="${fnamebase}.solid_kmers_binary";
filePacked="${fnamebase}.solid_kmers_binary.packed";

cmd="$CosmoPack $fileKmer"
echo $cmd
$cmd

# cmd="cat $fileFasta $fileFasta2" 
# echo $cmd 
# $cmd > merged.fasta

# cmd="$dsk merged.fasta $k"
# echo $cmd
# $cmd

# fileKmer2="merged.solid_kmers_binary";
# filePacked2="merged.solid_kmers_binary.packed";

# cmd="$CosmoPack $fileKmer2"
# echo $cmd
# $cmd

cmd="$CosmoAddVerify $filePacked $fileFasta"
echo $cmd
$cmd

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

