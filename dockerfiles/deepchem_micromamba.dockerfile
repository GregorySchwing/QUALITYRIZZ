# Issue build command from main directory 
FROM mambaorg/micromamba:1.4.9 
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/deepchem.yml env.yml 
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean --all --yes 
RUN micromamba install -y -n base -c conda-forge procps-ng && micromamba clean --all --yes 
RUN micromamba shell init -s bash -p ~/micromamba

#ENV BASH_ENV=/root/.bashrc

# Set the default command to run when the container starts
CMD ["/bin/bash"]
