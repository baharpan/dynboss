
dsk="../dsk-1.6906/dsk"
CosmoPack="../cosmo-pack"
cosmodel="../cosmo-delete-add"
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

cmd="$cosmodel $filePacked"
echo $cmd
$cmd



