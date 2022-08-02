FROM continuumio/miniconda3

SHELL ["/bin/bash", "-c"]

ADD prolonged_infection.yaml /tmp/prolonged_infection.yaml

RUN conda env create -f /tmp/prolonged_infection.yaml -n wisc-prolonged-infection python=3.8 \
	&& conda clean --all -f --yes \
	&& conda init bash \
	&& echo "conda activate wisc-prolonged-infection" > ~/.bashrc

ENV PATH /opt/conda/envs/wisc-prolonged-infection/bin:$PATH

ENV CONDA_DEFAULT_ENV wisc-prolonged-infection

CMD ["/bin/bash"]
