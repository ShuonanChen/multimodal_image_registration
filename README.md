# Multimodal fluorescence microscopy image registration using morphological and spatial information

This repository contains all the code used to run our code, to register *in vivo* and *ex vivo* images. Note that we assume that the cells in the images are segmented properly (we used [Cellpose](https://cellpose.readthedocs.io/en/latest/command.html) with some training data which works great). 

## The files needed
1. *ex vivo* image (as the target)
2. *in vivo* image (as the source)
1. *ex vivo* segmented result
2. *in vivo* segmented result.

## The basic steps in our pipelines
1. Segmentation (not in this repo)
2. Find the potentially matching cells using ICP to incorporate the information on cell morphology as well as the context. 
3. Find the consensus set by finding consensus matches and the maximum-rotation-set. 
4. Include more cells to learn the scale of the transformation by neighbor mathching. 
5. Fine-tuning using iterative non-rigid transformation based on the phase correlation

## The steps to run the code
All our functions are available from `./functions`. Simply run the jupyter notebook `./run_registration.ipynb`. This notebook runs the above steps and generate the results of 

## Version information
1. Python 3.7
2. MATLAB2022b


