FROM us.gcr.io/broad-dsde-methods/terra-jupyter-minimal-gpu-base:0.0.2
# Build off of image https://github.com/sjfleming/terra-docker/blob/sf_minimal_base/terra-jupyter-minimal-gpu-base/CHANGELOG.md
# terra-jupyter-minimal-gpu-base image
USER root

ENV PATH=$PATH:/app

COPY requirements.txt .
RUN pip3 install --break-system-packages -r requirements.txt

WORKDIR /app