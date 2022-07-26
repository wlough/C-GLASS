FROM dasteck/cglass:latest
WORKDIR /build
RUN git clone --recursive --single-branch --branch Daniel_Temporary https://github.com/Betterton-Lab/C-GLASS.git . &&\
    ./install.sh -otI
RUN rm -rf /build
WORKDIR /mnt
ENTRYPOINT []
