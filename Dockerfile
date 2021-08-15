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
RUN apt-get -y install r-base
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('qvalue')"


COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN mkdir -p /project /nl /mnt /share
ENV PATH /opt/conda/envs/dolphinnext-PSI-Sigma-4.0/bin:$PATH

# Install PSI-Sigma
RUN wget https://github.com/wososa/PSI-Sigma/archive/v1.9p.tar.gz && \
    tar -xzf v1.9p.tar.gz && mv PSI-Sigma-1.9p /usr/local/bin/PSI-Sigma-1.9p
ENV PATH /usr/local/bin/PSI-Sigma-1.9p:$PATH

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
RUN apt-get install -y git make g++ gcc python wget libgsl0-dev

# 3. Install PDL::GSL	
RUN cpanm PDL::LiteF
RUN cpanm PDL::GSL::CDF	
RUN cpanm PDL::Stats	
RUN cpanm Statistics::Multtest	
RUN cpanm Statistics::R	

# Install MEME
RUN apt-get install -y autoconf automake libtool ghostscript
RUN mkdir /opt/meme
ADD http://meme-suite.org/meme-software/5.3.3/meme-5.3.3.tar.gz /opt/meme
WORKDIR /opt/meme/
RUN tar zxvf meme-5.3.3.tar.gz && rm -fv meme-5.3.3.tar.gz
RUN cd /opt/meme/meme-5.3.3 && \
        ./configure --prefix=/opt  --enable-build-libxml2 --enable-build-libxslt && \
        make && \
        make install && \
        rm -rfv /opt/meme
ENV PATH="/opt/bin:${PATH}"

