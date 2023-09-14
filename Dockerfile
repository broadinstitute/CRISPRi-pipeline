FROM us.gcr.io/landerlab-atacseq-200218/gcloud-samtools:0.1
# TODO: change image
WORKDIR /app
ENV PATH=$PATH:/app

COPY requirements.txt .
RUN pip3 install --break-system-packages -r requirements.txt

COPY CRISPRi_pipeline ./CRISPRi_pipeline