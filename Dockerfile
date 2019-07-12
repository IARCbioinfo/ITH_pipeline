################## BASE IMAGE ######################
FROM nfcore/base

################## METADATA ######################

LABEL base_image="nfcore/base"
LABEL version="1.0"
LABEL software="svaba-nf"
LABEL software.version="1.0"
LABEL about.summary="Container image containing all requirements for ITH_pipeline"
LABEL about.home="http://github.com/delhommet/ITH_pipeline"
LABEL about.documentation="http://github.com/delhommet/ITH_pipeline/README.md"
LABEL about.license_file="http://github.com/delhommet/ITH_pipeline/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER Tiffany Delhomme <delhommet@students.iarc.fr>

################## INSTALLATION ######################
RUN conda install -c bioconda bcftools=1.7
RUN conda install -c bioconda samtools=1.7

# installation of python libraries
RUN pip install multiprocess && pip install subprocess && pip install pandas && pip install seaborn

# installation of bnpy
RUN cd ~ && git clone https://michaelchughes@bitbucket.org/michaelchughes/bnpy-dev/ && export BNPYROOT=~/bnpy-dev
RUN pip install numpy && pip install matplotlib

# environment variable to use correctly bnpy
ENV BNPYROOT=~/bnpy-dev
ENV PYTHONPATH=${PYTHONPATH}:~/bnpy-dev
ENV BNPYOUTDIR=~/nbpy-dev-results

# install gurobi
RUN cd ~ && wget https://packages.gurobi.com/8.1/gurobi8.1.1_linux64.tar.gz && tar -zxvf gurobi8.1.1_linux64.tar.gz
ENV PATH=${PATH}:~/gurobi811/linux64/bin
ENV GRB_LICENSE_FILE=~/gurobi811/gurobi.lic
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:~/gurobi811/linux64/lib"

# install hatchet
RUN cd ~ && git clone https://github.com/raphael-group/hatchet
RUN cd hatchet/ && mkdir build && cd build/ && cmake .. && make
