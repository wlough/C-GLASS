FROM shfiorenza/cglass_base:latest
WORKDIR /build
RUN git clone --recursive --single-branch --branch nab_xfer https://github.com/Betterton-Lab/C-GLASS.git . &&\
    ./install.sh -I
RUN rm -rf /build
WORKDIR /mnt
ENTRYPOINT []
