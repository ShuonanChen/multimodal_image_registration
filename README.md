# Multimodal fluorescence microscopy image registration using morphological and spatial information

This repository contains all the code used to run our code, to register *in vivo* and *ex vivo* images. Note that we assume that the cells in the images are segmented properly (we used [Cellpose](https://cellpose.readthedocs.io/en/latest/command.html) with some training data which works great). 

## The files needed (assuming pickle files)
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
All our functions are available from `./functions`. 
Simply run the jupyter notebook `./run_registration.ipynb`. This notebook runs the above steps and output the following as the results:
1. transformation parameters learned
2. registered image 
3. number of the matching cells. 


## Version information and dependencies
1. Python 3.9.7
2. MATLAB2021b
3. Python package dependencies:

```
# platform: linux-64

scipy==1.7.3
numpy==1.21.2
scikit-image==0.18.3
matplotlib==3.5.0
connected-components-3d==3.8.0
open3d==0.14.1
opencv-python==4.5.5.62
pandas==1.3.5
scikit-image==0.18.3
tqdm==4.62.3
```

