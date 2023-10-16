#! /bin/bash
#
# Run to launch a docker container named "cglass_latest".
# Provide the flag '-b' to force rebuilding the project image:
# >> ./run_docker.sh -b
# Provide the flag '-x' to launch/build experimental image from jmm/experimental branch
# >> ./run_docker.sh -bx

show_help() {
    echo "Without additional options, $0 launches a docker container named cglass_latest to run in the background"
    echo "USAGE:"
    echo "  $0 [-hbx]"
    echo "OPTIONS:"
    echo "  -h      show this menu"
    echo "  -b      force rebuild of C-GLASS Docker image"
    echo "  -x      build/launch experimental version of C-GLASS Docker image"
}

# Reset in case getopts has been used previously in the shell.
OPTIND=1         
experimental=false
build=false
while getopts "h?bx" opt; do
    case "$opt" in
    h|\?)
        show_help
        exit 0
        ;;
    b)  
        build=true
        ;;
    x)  
        experimental=true
        ;;
    esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift

if $build; then
    if $experimental; then
        echo "Building experimental docker image"
        docker build --no-cache -f docker/Dockerfile_experimental -t jeffmm/cglass:experimental docker
    else
        echo "Building docker image"
        docker build --no-cache -t shfiorenza/cglass:latest docker
    fi
elif $experimental; then
    echo "Launching experimental C-GLASS docker container"
    docker run --rm -itd -v "${PWD}":/mnt --name "cglass_experimental" jeffmm/cglass:experimental bash
# Otherwise, just start up the containers
else
  echo "Launching C-GLASS docker container"
  docker run --rm -itd -v "${PWD}":/mnt --name "cglass" shfiorenza/cglass:latest bash
fi
