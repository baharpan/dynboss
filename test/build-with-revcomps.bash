dsk="../dsk-1.6906/dsk"
CosmoPack="../cosmo-pack"
cosmodel="../cosmo-build-dyn"

if [ $# -lt 2 ]; then
    echo "Usage: $0 <fastafilebase> <k>"
    exit
fi

k=$2

cd ..
rm cosmo-pack
make cosmo-pack revcomps=1
rm cosmo-build-dyn
make cosmo-build-dyn revcomps=1
cd test

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

cmd="$cosmodel $filePacked"
echo $cmd
$cmd



