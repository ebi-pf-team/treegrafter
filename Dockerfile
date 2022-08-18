FROM ubuntu:latest

ARG branch=main

RUN apt update
RUN apt install -y --no-install-recommends bison build-essential ca-certificates cmake flex python3 unzip wget zlib1g-dev

RUN wget -O /opt/hmmer-3.3.1.tar.gz http://eddylab.org/software/hmmer/hmmer-3.3.1.tar.gz
RUN tar -C /opt/ -zxf /opt/hmmer-3.3.1.tar.gz
WORKDIR /opt/hmmer-3.3.1
RUN ./configure
RUN make
RUN make install

RUN wget -O /opt/epa-ng-0.3.7-zip https://github.com/Pbdas/epa-ng/archive/refs/tags/v0.3.7.zip
RUN unzip -d /opt/ /opt/epa-ng-0.3.7-zip
WORKDIR /opt/epa-ng-0.3.7
RUN make
RUN cp bin/epa-ng /usr/local/bin/

RUN wget -O /opt/treegrafter.zip https://github.com/ebi-pf-team/treegrafter/archive/refs/heads/${branch}.zip
RUN unzip -d /opt/ /opt/treegrafter.zip
RUN mv /opt/treegrafter-* /opt/treegrafter

WORKDIR /
RUN rm -rf /opt/epa-ng-0.3.7* /opt/hmmer-3.3.1* /opt/treegrafter.zip

ENTRYPOINT ["/opt/treegrafter/treegrafter.sh"]
