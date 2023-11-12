# Issue build command from main directory 
FROM mambaorg/micromamba:1.4.9 
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/openff-toolkit.yml /tmp/env.yml 
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean --all --yes 
RUN micromamba install -y -n base -c conda-forge procps-ng && micromamba clean --all --yes 

RUN mamba init

ENV BASH_ENV=/root/.bashrc

# Set the default command to run when the container starts
CMD ["/bin/bash"]
