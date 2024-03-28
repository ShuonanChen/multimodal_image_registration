import time, sys
import numpy as np
import pandas as pd
import scipy.ndimage as ndi
from skimage import color, morphology, measure
from skimage.transform import downscale_local_mean, resize, rescale, warp
from skimage.registration import phase_cross_correlation
from scipy import sparse
import matplotlib.pyplot as plt
from cellregister import ncc
import tiling

def learn_and_apply_deformable(im1_happy,im2,vec_ds, patch_sizes = (100,100,100), patch_boarder=(20,20,20), deform_kernel = 100):
    time1=time.time()
    im1_foo = downscale_local_mean(im1_happy,(vec_ds,vec_ds,vec_ds))  # to be transforme (zstack)
    im2_foo = downscale_local_mean(im2,(vec_ds,vec_ds,vec_ds)) # 
    ncc_1 = np.around(ncc(im1_foo,im2_foo), 3)
    tile_szs = (patch_sizes[0]//vec_ds,patch_sizes[1]//vec_ds,patch_sizes[2]//vec_ds)
    border_szs = (patch_boarder[0]//vec_ds,patch_boarder[1]//vec_ds,patch_boarder[2]//vec_ds)
    tiles = tiling.tile_up_nd(im1_foo.shape, tile_szs, border_szs=border_szs)
    ksz=deform_kernel//vec_ds
    out = im1_foo.copy()
    ncc_list = []
    vec_field_smooth_list = []  # so that we can apply to the binary image later!
    for i in range(30): 
        if i%5==0:
            print(f'iteration {i}...')
        vec_field = np.zeros(out.shape + (3,))
        pix_counts = np.zeros_like(out)
        for t in tiles:        
            im1_t_f = im2_foo[t.look.as_slices]
            im2_t_f = out[t.look.as_slices]    
            try:
                shift,_,_ = phase_cross_correlation(im1_t_f, im2_t_f)    
            except:
                shift = np.zeros(3)
            if np.abs(shift).max() > 30//vec_ds:
                vec_field[t.put.as_slices] += np.zeros(3)[None,None,None,:]
            else:
                shifted = ndi.shift(im2_t_f, shift)
                delta_r = ncc(im1_t_f,shifted) - ncc(im1_t_f,im2_t_f)
                if delta_r>0:  
                    vec_field[t.look.as_slices] += delta_r*shift[None,None,None,:]
                    pix_counts[t.look.as_slices] += delta_r
                else:
                    vec_field[t.put.as_slices] += np.zeros(3)[None,None,None,:]
        pix_counts[pix_counts==0] = 1e-5
        vec_field = vec_field/pix_counts[:,:,:,None]     
        vec_field[np.isnan(vec_field)] = 0    
        vec_field_smooth = vec_field.copy()
        for coord in range(3):
            vec_field_smooth[:,:,:,coord] = ndi.gaussian_filter(vec_field[:,:,:,coord], ksz)
        vec_field_smooth_list.append(vec_field_smooth)

        shape = out.shape
        mapz_base, mapx_base, mapy_base = np.meshgrid(np.arange(shape[0]),np.arange(shape[1]),np.arange(shape[2]), indexing='ij')
        mapz=mapz_base-vec_field_smooth[:,:,:,0]
        mapx=mapx_base-vec_field_smooth[:,:,:,1]
        mapy=mapy_base-vec_field_smooth[:,:,:,2]

        img_de = warp(out,np.array((mapz,mapx,mapy)), order = 3)
        ncc_list.append(ncc(im2_foo,img_de))
        if i>3 and ncc_list[-1]<=np.mean(np.array(ncc_list[:3])):
            print('stopping since NCC is not moving up')
            break;   
        if i>5 and ncc_list[-1]==ncc_list[-3]:
            print('stopping since NCC is not moving up')
            break;
        out = img_de        
    img_de = out  
    im1_de = rescale(img_de,(vec_ds,vec_ds,vec_ds), order = 3)
    
    plt.figure(figsize=(4,3))
    plt.plot(ncc_list, '-o')
    plt.axhline(ncc_1, color = 'red')
    plt.ylabel('NCC')
    plt.xlabel('iteration')
    out = [img_de,im1_de,ncc_list,vec_field_smooth_list]
    return(out)
