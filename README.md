
# ProkGS+ : Prokaryotic Gene Structure plus

A bioinformatics tool that can help users to cross validate functional annotations of structurally characterized proteins from the Protein Data Bank (PDB) by mapping them to their respective genomes.
 


## Introduction
ProkGS is a simple,easy-to-use tool that helps user validate if the functional annotation of a protein is consistent when mapped to its respective genome. It does so by following a meticulous pipeline that maps the PDB protein to its source genome and obtains the GFF(general feature format) file which contains the CDS(coding sequence) features from which the annotaion of the protein function can be extracted.
## Requirements
To run ProkGS+ the only requirement is to install Docker on your system. Click on the following links to download Docker CE engine for your system.

NOTE: ProkGS is currently available for the following distributions. To use it on Windows or MacOs please install Docker Desktop: [Click here](https://docs.docker.com/desktop/) 

 - [Ubuntu](https://docs.docker.com/engine/install/ubuntu/)
 - [Debian](https://docs.docker.com/engine/install/debian/)
 - [RHEL](https://docs.docker.com/engine/install/rhel/)
 - [Fedora](https://docs.docker.com/engine/install/fedora/)
 - [CentOs](https://docs.docker.com/engine/install/centos/)

If the user wants to use the GEMINI feature of the pipeline, then they need an API key.

### How to generate an API key for free:
Follow this 2minute tutorial:
-[Click here](https://youtu.be/c3GWAbxyr3A?si=ZmG0nyG62PT0BEHY)

                              OR 

Create or login to your google account and click on "get API key" [here](https://aistudio.google.com/app/apikey)

## How to Use
To start running ProkGS+ after docker installation, create a folder named "DOCKER" on your Desktop, and copy the chiranth.py,ggplot.R and Dockerfile onto this directory/folder. Once this is done follow the Steps outlined below: 

### Step 1:
Pull the docker image
```bash
  docker pull chiranth27/prokgs_plus:latest
```
### Step 2:
Check if the image exists
```bash
  docker images
```
There are two ways a user can use ProkGS+:

#### Without API KEY:
Run this command if you have an API key:
```bash
  docker run -it --rm -v "$(pwd)/Results":/home/proKGS/Results chiranth27/prokgs_plus

```
#### With API KEY:
```bash
 docker run -it --rm -e GEMINI_API_KEY=ABC123 -v "$(pwd)/Results":/home/proKGS/Results chiranth27/prokgs_plus

```
Replace ABC123 with your Gemini API key.

Both the commands prompts docker to run in an interactive mode and opens up a terminal for the user where they can enter their 4 letter PDB ID


## Outputs

ProkGS+ outputs 4 files:
- Stream.txt- which contains all the information needed to plot the final graph
- Gff file- file from which the gene annoation was obtained.
- PNG image- An image that gives gene neighbourhood context for better understanding the gene in question.[Click here for an example image](https://ibb.co/KcT4717X)
- csv file- file that gives information about the pdb id, the genome accession, gene and pdb annoation and a flag(if run using an API key)
- **Disclaimer**:Generative AI is experimental and any sound decisions or conclusions should not be made from its classification which is a  part of the ProkGS+ pipeline.
## Authors
- Chiranth M - [mail](chiranth.m2001@gmail.com)
- Dr. S Thiyagarajan -[mail](sthiyaga@ibab.ac.in)


