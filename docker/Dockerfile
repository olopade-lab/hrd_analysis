FROM r-base:4.0.5

RUN apt-get update \
	&& apt-get install -y apt-utils libcurl4-openssl-dev libxml2-dev libfontconfig1-dev libssl-dev

RUN R -e "install.packages(c('optparse', 'BiocManager', 'randomForest', 'devtools'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "BiocManager::install(c('BSgenome', 'BSgenome.Hsapiens.UCSC.hg19'))"

RUN R -e "BiocManager::install(c('BSgenome', 'BSgenome.Hsapiens.UCSC.hg38'))"

RUN R -e "devtools::install_github('https://github.com/UMCUGenetics/mutSigExtractor/')"

RUN R -e "devtools::install_github('https://github.com/UMCUGenetics/CHORD/')"
