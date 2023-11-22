FROM ubuntu:18.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites and R
RUN apt update && \
    ln -fs /usr/share/zoneinfo/Europe/Dublin /etc/localtime && \
    apt-get install -y build-essential wget unzip r-base libssl-dev libxml2-dev libcurl4-openssl-dev

# Install devtools
RUN R -e "install.packages('curl', repos = 'http://cran.us.r-project.org'); library(curl)"
RUN R -e "install.packages('httr', repos = 'http://cran.us.r-project.org'); library(httr)"
RUN R -e "install.packages('usethis', repos = 'http://cran.us.r-project.org'); library(usethis)"
RUN R -e "install.packages('devtools', repos = 'http://cran.us.r-project.org'); library(devtools)"
RUN R -e "install.packages('tidyverse', repos = 'http://cran.us.r-project.org'); library(tidyverse)"
RUN R -e "install.packages('vroom', repos = 'http://cran.us.r-project.org'); library(vroom)"

# Install corncob
RUN R -e "library(devtools); devtools::install_github('statdivlab/corncob')"
