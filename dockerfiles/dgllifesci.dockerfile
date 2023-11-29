# Issue build command from main directory 
FROM nvcr.io/nvidia/dgl:23.09-py3
RUN apt-get -y update && pip install rdkit
RUN apt-get -y update && pip install dgllife
# Set the default command to run when the container starts
CMD ["/bin/bash"]
