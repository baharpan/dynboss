cosmodel="../bin/cosmo-delete-verify"
k=31

if [ $# -lt 2 ]; then
    echo "Usage: $0 <.dbg file> <.fasta file>"
    exit
fi

graph="$1"
fastaFile="$2"

cmd="$cosmodel $graph $fastaFile"
echo $cmd
$cmd
