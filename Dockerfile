FROM ubuntu:latest

ARG DEBIAN_FRONTEND=noninteractive

# RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-cran-randomforest python3.10 python3-pip python3-setuptools python3-dev
RUN apt-get update && apt-get install -y --no-install-recommends \
 build-essential \
 r-base \
 r-base-dev \
 python3.10 \
 python3-pip \
 python3-setuptools \
 python3-dev

RUN apt-get install -y --no-install-recommends \
 libxml2-dev \
 libpq-dev \
 libssl-dev \
 libcurl4-openssl-dev \
 libfontconfig1-dev \
 libharfbuzz-dev \
 libfribidi-dev \
 libgit2-dev \
 libfreetype6-dev \
 libpng-dev \
 libtiff5-dev \
 libjpeg-dev \
 libcairo2-dev \
 libxt-dev \
 libproj-dev \
 bowtie2 \
 subread

 RUN apt-get -y install r-cran-devtools libgit2-dev

# Copy Necessary Python and R Requirement Files and Install
COPY /home/requirements.txt /home/requirements.txt
RUN pip3 install -r /home/requirements.txt
COPY /home/requirements.r /home/requirements.r
RUN Rscript /home/requirements.r

# Install SAMTOOLS
COPY /home/installations/samtools-1.12.tar.bz2 /home/installations/samtools-1.12.tar.bz2
WORKDIR /home/installations/
RUN tar xvjf /home/installations/samtools-1.12.tar.bz2
WORKDIR /home/installations/samtools-1.12/
RUN ./configure
RUN make
RUN make install
WORKDIR /home

# Set up the User Interface
RUN apt-get -y install sudo

RUN useradd --create-home -s /bin/bash cs181-user && \
  echo "cs181-user ALL=(ALL:ALL) NOPASSWD: ALL" > /etc/sudoers.d/cs181-init

# RUN mkdir /home/cs181-user
# RUN chown cs181-user /home/cs181-user

USER cs181-user
RUN rm -f ~/.bash_logout

WORKDIR /home

RUN sudo apt-get -y install git

CMD ["/bin/bash", "-l"]