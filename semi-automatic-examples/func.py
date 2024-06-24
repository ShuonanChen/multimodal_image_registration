import numpy as np
import scipy.ndimage as ndi


def trans(source,T_dict,o_shape):  # original when source image is of size (c,z,x,y)
    assert(len(source.shape)==4)  # this is (c,z,x,y) where c==2 for now. z can be 8,16,24....
    C = source.shape[0]
    B = np.eye(4)
    for l in T_dict:
        for k,v in l.items():
            if k=='bhat':
                B = B@((np.c_[v, np.array((0,0,0,1))]))
            if k=='scale':
                B[:,:3] *= v  
    R_3 = (np.linalg.inv(B[:3,:3])).T
    offset_3 = -B[-1,:-1]@np.linalg.inv(B[:3,:3])
    print('running rigid..')
    transformed_all = np.array([ndi.affine_transform(source[c].astype('float32'), R_3, offset = offset_3,
                                    output_shape = o_shape, order=1) for c in range(C)])
    return(transformed_all)


def trans_w_deform(source,T_dict,o_shape):  # original when source image is of size (c,z,x,y)
    '''o_shape = (source.shape[1],) + exvivos[section_id][0].shape[1:] -- # (z,x,y)'''
    assert(len(source.shape)==4)  # this is (c,z,x,y) where c==2 for now. z can be 8,16,24....
    all_vec_f_3 = np.zeros(tuple(o_shape) + (3,))
    C = source.shape[0]
    B = np.eye(4)
    for l in T_dict:
        for k,v in l.items():
            if k=='bhat':
                B = B@((np.c_[v, np.array((0,0,0,1))]))
            if k=='scale':
                B[:,:3] *= v  
            if k=='vec_field_total':
                vf = np.stack([pad3d(v[...,c],out_shape = o_shape) for c in range(3)], -1)
                all_vec_f_3 += vf
    R_3 = (np.linalg.inv(B[:3,:3])).T
    offset_3 = -B[-1,:-1]@np.linalg.inv(B[:3,:3])
    print('running rigid..')
    transformed_all = np.array([ndi.affine_transform(source[c].astype('float32'), R_3, offset = offset_3,
                                    output_shape = o_shape, order=1) for c in range(C)])  
    mapz_base, mapx_base, mapy_base = np.meshgrid(np.arange(o_shape[0]),np.arange(o_shape[1]), np.arange(o_shape[2]),indexing='ij')
    mapz=mapz_base-all_vec_f_3[:,:,:,0]
    mapx=mapx_base-all_vec_f_3[:,:,:,1]
    mapy=mapy_base-all_vec_f_3[:,:,:,2]
    print('running deformable..')
    deformed_all = np.array([warp(transformed_all[c],np.array((mapz,mapx,mapy)), order = 1) for c in range(C)])
    return(deformed_all)






def pad3d(a,out_shape):
    '''only deal with the case out_shape is smaller in all three axis!'''
    assert(len(a.shape)==3)
    assert(len(a.shape)==len(out_shape))    
    out =a[:out_shape[0],:out_shape[1],:out_shape[2]]
    return out