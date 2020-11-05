#!/bin/bash
set -eo pipefail

# adjust these parameters to your needs
DOCKER_IMAGE=ljelonek/bakta:dev
DEFAULT_DBPATH=$PWD/db

args=( "$@" )
argcount=${#args[@]}

# handle help parameter separately
if [[ $argcount -eq 0 || " ${args[@]} " =~ " --help " || " ${args[@]} " =~ " -h " ]]; then
    docker run -it --rm \
               --user $(id -u):$(id -g) \
               $DOCKER_IMAGE
    exit $?
fi


# run bakta docker
# decisions:
# all parameters are passed through, except genome, db and output
# they are replaced with their absolute paths and are mounted to fixed
# positions in the container
#
# genome -> /bakta/genome.fas
# output -> /bakta/output
# db -> /bakta/db

DB=
OUTPUT=
GENOME=$(realpath ${args[$argcount-1]})
args[$((argount-1))]=/bakta/genome.fas

for i in `seq 0 $((argcount-1))`; do
    VAL=${args[i]}
    j=$((i+1))
    if [[ "$VAL" == "--db" || "$VAL" == "-d" ]]; then
        DB=${args[j]}
        DB=$(realpath $DB)
        args[j]=/bakta/db
    fi
    if [[ "$VAL" == "--output" || "$VAL" == "-o" ]]; then
        OUTPUT=${args[j]}
        OUTPUT=$(realpath $OUTPUT)
        args[j]=/bakta/output
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
echo "* bakta docker wrapper       *"
echo "******************************"
echo "*       # Arguments: " $argcount
echo "* Database location: " $DB
echo "*   Genome location: " $GENOME
echo "*  Output location : " $OUTPUT
echo "******************************"
CMD=$(cat <<-END
    docker run -it --rm \
    --user $(id -u):$(id -g)
    -v $DB:/bakta/db:ro \
    -v $OUTPUT:/bakta/output:rw \
    -v $GENOME:/bakta/genome.fas:ro \
    $DOCKER_IMAGE ${args[@]}
END
)
echo "* Commandline: " $CMD
echo "******************************"
$CMD
