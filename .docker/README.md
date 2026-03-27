# Dockerfile

This directory contains a Dockerfile for creating a Docker image. 

This Docker image contains all Conda environments for each rule, i.e. the whole workflow is run in one image.

These images are shared via [Docker Hub](https://hub.docker.com/repository/docker/niekwit/smallrna-seq/general) and are generated as follows (from directory with workflow code):

```shell
$ snakemake --containerize > Dockerfile
$ sudo docker build -t niekwit/smallrna-seq:{VERSION} .
$ sudo docker login
$ sudo docker push niekwit/smallrna-seq:{VERSION}
```

When `Snakemake` is run with `--use-apptainer True`, the workflow will automatically pull the latest image from Docker Hub and convert it to an Apptainer image on the fly.