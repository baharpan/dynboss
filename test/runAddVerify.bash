CosmoAddVerify="../cosmo-add-verify"
k=31


if [ $# -lt 2 ]; then
    echo "Usage: $0 <.dbg file> <.fasta file that DynamicBOSS is built on>"
    exit
fi

graph="$1"
fastaFile="$2"

cmd="$CosmoAddVerify $graph $fastaFile"
echo $cmd
$cmd
