
dsk="../dsk-1.6906/dsk"
CosmoPack="../cosmo-pack"
CosmoAddRead="../cosmo-add-reads"
k=64

#OutFile = "./cmpResultsCosmoBowe.txt"

if [ $# -lt 1 ]; then
    echo "Usage: $0 <fastafilebase (graph)> <fastafilebase (reads to add)>"
    exit
fi

fnamebase="$1"
fileFasta="$1.fasta"

fnamebase2="$2"
fileFasta2="$2.fasta"

fileKmer="${fnamebase}.solid_kmers_binary";
if [ -f "$fileKmer" ]; then
   echo "DSK already ran."
else
    cmd="$dsk $fileFasta $k"
    echo $cmd
    $cmd
fi

filePacked="${fnamebase}.solid_kmers_binary.packed";

if [ -f "$filePacked" ]; then
   echo "Kmers already packed."
else
    cmd="$CosmoPack $fileKmer"
    echo $cmd
    $cmd
fi
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

cmd="$CosmoAddRead $filePacked $fileFasta2"
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

