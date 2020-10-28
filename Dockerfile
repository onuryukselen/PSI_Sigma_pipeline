FROM ubuntu:16.04
LABEL author="onur.yukselen@umassmed.edu"  description="Docker image containing all requirements for the PSI-Sigma pipeline"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc
    
# configure image 
RUN apt-get -y update 
RUN apt-get -y install software-properties-common build-essential
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran40/'
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
#RUN apt-key adv --keyserver pgp.mit.edu --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-get -y install apt-transport-https
RUN apt-get -y update

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN mkdir -p /project /nl /mnt /share
ENV PATH /opt/conda/envs/dolphinnext-PSI-Sigma-1.0/bin:$PATH

# Install PSI-Sigma
RUN wget https://github.com/wososa/PSI-Sigma/archive/v1.9j.tar.gz && \
    tar -xzf v1.9j.tar.gz && mv PSI-Sigma-1.9j /usr/local/bin/PSI-Sigma-1.9j
ENV PATH /usr/local/bin/PSI-Sigma-1.9j:$PATH

# Install compiler and perl stuff
RUN apt-get install --yes \
 build-essential \
 gcc-multilib \
 apt-utils \
 perl \
 expat \
 libexpat-dev 

# SET PERL5LIB
ENV PERL5LIB="/usr/local/lib/x86_64-linux-gnu/perl/5.22"
# 1. Install cpanm
RUN cpan App::cpanminus
# RUN cpanm PDL::LiteF
# RUN cpanm PDL::Stats

# 2. Install GSL 
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y git make g++ gcc python wget


ENV GSL_TAR="gsl-2.4.tar.gz"
ENV GSL_DL="ftp://ftp.gnu.org/gnu/gsl/$GSL_TAR"

ENV GSL_ROOT="/usr/local/lib/x86_64-linux-gnu/perl/5.22"
ENV LD_LIBRARY_PATH="$GSL_ROOT/lib:$LD_LIBRARY_PATH"

ENV GSL_INC="/usr/local/lib/x86_64-linux-gnu/perl/5.22/include"

RUN wget -q $GSL_DL \
    && tar zxvf $GSL_TAR \
    && rm -f $GSL_TAR \
    && cd gsl-2.4 \
    && ./configure --prefix=/usr/local/lib/x86_64-linux-gnu/perl/5.22 \
    && make -j 4 \
    && make install

# 3. Install PDL::GSL
RUN cpanm PDL::LiteF
RUN cpanm PDL::Stats

RUN wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz \
    && tar zxvf gsl-2.4.tar.gz \
    && cd gsl-2.4 \
    && ./configure \
    && make \
    && make install 

RUN cpanm PDL::GSL::CDF
RUN cpanm Statistics::Multtest