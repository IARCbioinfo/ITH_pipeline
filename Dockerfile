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

RUN apt-get update -y && \
	DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y build-essential # to have gcc compiler for seaborn

RUN cd ~ && wget http://www.cmake.org/files/v3.2/cmake-3.2.2.tar.gz && tar xf cmake-3.2.2.tar.gz && cd cmake-3.2.2 && ./configure && make
RUN chmod +x ~/cmake-3.2.2/bin/cmake && ln -s ~/cmake-3.2.2/bin/cmake /usr/bin

################## HATCHET ###########################
RUN conda install -c bioconda bcftools=1.7
RUN conda install -c bioconda samtools=1.7

# installation of python libraries
RUN pip install multiprocess && pip install pandas && pip install seaborn && pip install matplotlib
RUN pip install matplotlib-venn #in order to run the BBeval script

# installation of bnpy
RUN cd ~ && git clone https://michaelchughes@bitbucket.org/michaelchughes/bnpy-dev/

# environment variable to use correctly bnpy
ENV BNPYROOT=~/bnpy-dev
ENV PYTHONPATH=${PYTHONPATH}:~/bnpy-dev
ENV BNPYOUTDIR=~/nbpy-dev-results

# install gurobi
RUN cd ~ && wget https://packages.gurobi.com/8.1/gurobi8.1.1_linux64.tar.gz && tar -zxvf gurobi8.1.1_linux64.tar.gz
RUN cd ~/gurobi811/linux64/src/build/ && make && cp libgurobi_c++.a ../../lib/
ENV PATH=${PATH}:~/gurobi811/linux64/bin
# next should be modified to get the license path in input
ENV GRB_LICENSE_FILE=~/gurobi811/gurobi.lic
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:~/gurobi811/linux64/lib"
ENV GUROBI_LIB=~/gurobi811/linux64/lib/libgurobi81.so

# install hatchet
# here should 1. modify FindGUROBI.cmake to add gurobi path 2. modify CMakeLists.txt by adding -pthread
RUN cd ~ && git clone https://github.com/raphael-group/hatchet
RUN cd ~/hatchet && sed -i '5s/""/"~\/gurobi811\/"/g' FindGUROBI.cmake && sed -i '7s/-std=c++11/-std=c++11 -pthread/g' CMakeLists.txt # second sed only in this dockerfile
RUN mkdir ~/hatchet/build && cd ~/hatchet/build/ && cmake .. \
        -DGUROBI_CPP_LIB=/root/gurobi811/linux64/lib/libgurobi_c++.a \
        -DGUROBI_INCLUDE_DIR=/root/gurobi811/linux64/include/ \
        -DGUROBI_LIB=/root/gurobi811/linux64/lib/libgurobi81.so && make
RUN sed -i '451,452s/.*/#&/' ~/hatchet/utils/ArgParsing.py # see issue: https://github.com/raphael-group/hatchet/issues/1

################# DECIFER ##########################
# install boost library
#RUN cd ~ && wget https://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.tar.bz2
#RUN tar --bzip2 -xf boost_1_61_0.tar.bz2
#RUN ./bootstrap.sh --prefix=~/usr && ./b2
#ENV BOOST_ROOT=~/boost_1_61_0

# install lemon
#RUN cd ~ && wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz
#RUN tar -zxvf lemon-1.3.1.tar.gz
#RUN cd lemon-1.3.1/ && mkdir build && cd build && cmake .. && make

#install decifer
#RUN cd ~ && git clone https://github.com/raphael-group/decifer && cd decifer
#RUN mkdir build && cd build

#RUN cmake -DLIBLEMON_ROOT=~/lemon-1.3.1/lemon -DCPLEX=OFF \
#-DGUROBI_INCLUDE_DIR=~/gurobi811/linux64/include \
#-DGUROBI_CPP_LIB=~/gurobi811/linux64/lib/libgurobi_c++.a \
#-DGUROBI_LIB=~/gurobi811/linux64/lib/libgurobi81.so ..
