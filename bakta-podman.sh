#!/bin/bash
set -eo pipefail

# ToDo: set final namespace
IMAGE=quay.io/biocontainers/bakta:1.8.2--pyhdfd78af_0
DEFAULT_DBPATH=$BAKTA_DB

args=( "$@" )
argcount=${#args[@]}

# handle help parameter separately
if [[ $argcount -eq 0 || " ${args[@]} " =~ " --help " || " ${args[@]} " =~ " -h " ]]; then
    podman run -it --rm $IMAGE bakta --help
    exit $?
fi

# run bakta podman
# decisions:
# all parameters are passed through, except genome, db and output
# they are replaced with their absolute paths and are mounted to
# the container

CWD=$(realpath .)
DB=
OUTPUT=
GENOME=$(realpath ${args[$argcount-1]})
args[$((argcount-1))]=$GENOME

for i in `seq 0 $((argcount-1))`; do
    VAL=${args[i]}
    j=$((i+1))
    if [[ "$VAL" == "--db" || "$VAL" == "-d" ]]; then
        DB=${args[j]}
        DB=$(realpath $DB)
        args[j]=$DB
    fi
    if [[ "$VAL" == "--output" || "$VAL" == "-o" ]]; then
        OUTPUT=${args[j]}
        OUTPUT=$(realpath $OUTPUT)
        args[j]=$OUTPUT
    fi
done;

if [[ -z "$DB" ]]; then
    DB=$DEFAULT_DBPATH
fi
if [[ -z "$OUTPUT" ]]; then
    OUTPUT=$PWD
fi
if [[ ! -d "$OUTPUT" ]]; then
    mkdir -p "$OUTPUT"
fi

echo "******************************"
echo "* bakta podman wrapper       *"
echo "******************************"
echo "*       # Arguments: " $argcount
echo "* Database location: " $DB
echo "*   Genome location: " $GENOME
echo "*  Output location : " $OUTPUT
echo "******************************"
CMD=$(cat <<-END
    podman run -it --rm \
    -v $DB:$DB:ro \
    -v $OUTPUT:$OUTPUT:rw \
    -v $GENOME:$GENOME:ro \
    --workdir $CWD \
    $IMAGE bakta --force ${args[@]}
END
)
echo "* Commandline: " $CMD
echo "******************************"
$CMD
