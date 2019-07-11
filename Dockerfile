################## BASE IMAGE ######################
FROM nfcore/base

################## METADATA ######################

LABEL base_image="nfcore/base"
LABEL version="1.0"
LABEL software="svaba-nf"
LABEL software.version="1.0"
LABEL about.summary="Container image containing all requirements for svaba-nf"
LABEL about.home="http://github.com/IARCbioinfo/svaba-nf"
LABEL about.documentation="http://github.com/IARCbioinfo/svaba-nf/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/svaba-nf/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER Tiffany Delhomme <delhommet@students.iarc.fr>

################## INSTALLATION ######################
RUN conda install -c bioconda bcftools=1.7
RUN conda install -c bioconda samtools=1.7

# installation of python libraries
pip install multiprocess
pip install subprocess

# installation of bnpy
RUN cd ~ && git clone https://michaelchughes@bitbucket.org/michaelchughes/bnpy-dev/ && export BNPYROOT=~/bnpy-dev
RUN pip install numpy && pip install matplotlib

# installation of eigen
RUN wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2
RUN tar -vxjf 3.3.7.tar.bz2 && cd eigen*
