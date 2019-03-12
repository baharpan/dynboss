
dsk="../dsk-1.6906/dsk"
CosmoPack="../cosmo-pack"
cosmodel="../cosmo-delete-verify"
k=64

#OutFile = "./cmpResultsCosmoBowe.txt"

if [ $# -lt 2 ]; then
    echo "Usage: $0 <fastafilebase> <fastafile2base>"
    exit
fi

fnamebase="$1"
fileFasta="$1.fasta"

fnamebase2="$2"
fileFasta2="$2.fasta"

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

cmd="$cosmodel $filePacked $fileFasta2"
echo $cmd
$cmd



