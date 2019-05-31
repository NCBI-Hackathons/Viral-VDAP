FROM continuumio/miniconda3 #need version

ADD environment.yaml /tmp/environment.yaml

# create conda environment from yaml file
RUN conda env create --name viral-vdap --file /tmp/environment.yaml

# activate conda environment
RUN echo "source activate viral-vdap" > ~/.bashrc

# add programs to our path
ENV PATH /opt/conda/envs/env/bin:$PATH
