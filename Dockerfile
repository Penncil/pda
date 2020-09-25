FROM ubuntu:bionic
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
apt-utils \
ed \
less \
locales \
vim-tiny \
wget \
ca-certificates \
apt-transport-https \
gsfonts \
gnupg \
libcurl4-openssl-dev \
inetutils-ping \
curl \
libxml2-dev \
xz-utils
# Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/" > /etc/apt/sources.list.d/cran.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
ENV R_BASE_VERSION 4.0.2
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
                libssl-dev \
		r-base=${R_BASE_VERSION}* \
		r-base-dev=${R_BASE_VERSION}* \
		r-recommended=${R_BASE_VERSION}* \
        && echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site 
RUN R -e "install.packages(c('Rcpp', 'RcppArmadillo', 'httr', 'data.table','rvest'))"
RUN useradd pda -u 1000\
        && mkdir /home/pda
RUN addgroup pda staff
COPY build/pda_1.0.tar.gz /home/pda/pda
RUN chown -R pda:pda /home/pda
WORKDIR /home/pda
RUN R -e "install.packages('./pda',repos = NULL, type='source')"
CMD ["bash"]
