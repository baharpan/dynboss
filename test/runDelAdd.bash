
cosmodel="../bin/cosmo-delete-add"
k=31

if [ $# -lt 2 ]; then
    echo "Usage: $0 <.dbg file> <number of kmers to delete and add>"
    exit
fi

graph="$1"
numKmers="$2"

cmd="$cosmodel $graph $numKmers"
echo $cmd
$cmd
