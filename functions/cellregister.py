import os,sys, time, pickle
import numpy as np 
import pandas as pd
import matplotlib.pylab as plt
from scipy.io import savemat, loadmat
from skimage import color, measure
from skimage.morphology import white_tophat,disk
from skimage.transform import downscale_local_mean, resize, rescale, warp
from scipy import sparse
from scipy.stats import zscore
from scipy.spatial import distance
import cc3d
import open3d as o3d
import copy
import scipy.ndimage as ndi
from matplotlib.colors import ListedColormap
from matplotlib import cm
from pprint import pprint

    
def angles(alpha):
    '''
    alpha is angle with radian unit, returns the angle degree. 
    alpha can be numpy array or just a scalar.
    '''
    return(np.around(180*alpha/np.pi,3))

def normalize(img):
    foo = img.ravel()
    out = foo-foo.mean()
    out = out/np.std(foo)
    return(out)


def ncc(im1p,im2p):
    foo1 = normalize(im1p)
    foo2 = normalize(im2p)
    out = (foo1@foo2)/np.product(foo1.shape)
    return(out)



def crop(img_i,img_j):
    '''
    arbitrary cropping function, 
    !! NOTE this should be applied to the original images !!
    the two images are aligned at the center, 
    returns two images that are both smaller than the original. 
    '''
    cent_i  = ((np.array((img_i.shape))-1)/2).astype(int)
    cent_j  = ((np.array((img_j.shape))-1)/2).astype(int)
    new_sz = np.min(np.array((img_i.shape, img_j.shape)), axis = 0)
    sls = [slice(x,y) for (x,y) in  zip(cent_i-((new_sz-1)/2).astype(int),cent_i+(new_sz/2).astype(int)+1)]
    im1p = img_i[sls[0],sls[1],sls[2]]
    sls = [slice(x,y) for (x,y) in  zip(cent_j-((new_sz-1)/2).astype(int),cent_j+(new_sz/2).astype(int)+1)]
    im2p = img_j[sls[0],sls[1],sls[2]]
    assert(im1p.shape == im2p.shape)
    return(im1p,im2p)


def get_orig_2(orig_full, bbox, add_edge):
    '''
    get the target position (specified by bbox) from orig_full, 
    add edge too. 
    '''
    bb_with_edge = np.clip(bbox+np.c_[(-add_edge,add_edge)], 
                           np.zeros_like(bbox),
                           np.c_[orig_full.shape,orig_full.shape]).astype(int)
    cen = bbox[:,0] + (np.ptp(bbox, axis = 1)/2)
    new_radius = np.min(np.abs(bb_with_edge-cen[:,None]), axis = 1)
    bb_with_edge_fixed = (np.c_[cen-new_radius, cen+new_radius]).astype(int)

    sls = [slice(x,y) for (x,y) in bb_with_edge_fixed]
    orig_i = orig_full[sls[0],sls[1],sls[2]]
    return(orig_i)
    
    

def crop_w_given_sz(img_i, output_size):
    '''
    arbitrary cropping function, the image is cropped to the given output_size and returned
    !! NOTE this should be applied to the original images !! raises an error when the image size is smaller already    
    
    input:
        img_i: 3d size (m0,m1,m2)
    output:
        img_out: 3d size is equal to output_size. 
    '''
    assert(len(output_size)==3)
    assert( (np.array(img_i.shape)>=np.array(output_size)).all() )
    
    cent_i  = ((np.array((img_i.shape))-1)/2).astype(int)
    cent_j  = ((np.array((output_size))-1)/2).astype(int)
    new_sz = np.min(np.array((img_i.shape, output_size)), axis = 0)
    sls = [slice(x,y) for (x,y) in  zip(cent_i-((new_sz-1)/2).astype(int),cent_i+(new_sz/2).astype(int)+1)]
    im1p = img_i[sls[0],sls[1],sls[2]]
#     sls = [slice(x,y) for (x,y) in  zip(cent_j-((new_sz-1)/2).astype(int),cent_j+(new_sz/2).astype(int)+1)]
#     im2p = img_j[sls[0],sls[1],sls[2]]
    assert(im1p.shape == tuple(output_size))
    return(im1p)


# code sourse: https://learnopencv.com/rotation-matrix-to-euler-angles/ --  Checks if a matrix is a valid rotation matrix.
def isRotationMatrix(R):
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype = R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6

# Calculates rotation matrix to euler angles -- The result is the same as MATLAB!!
def rotationMatrixToEulerAngles(R):
    import math
    assert(isRotationMatrix(R))
    sy = math.sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])
    singular = sy < 1e-6
    if  not singular :
        x = math.atan2(R[2,1] , R[2,2])
        y = math.atan2(-R[2,0], sy)
        z = math.atan2(R[1,0], R[0,0])
    else :
        x = math.atan2(-R[1,2], R[1,1])
        y = math.atan2(-R[2,0], sy)
        z = 0
    return np.array([z, y, x])


def generate_stuff(im2,im2_bi,add_edge,ds_sz,outputshape,ds_sz_aff):
    regionprops_f = measure.regionprops(im2_bi)
    bboxes = np.array([np.array(prop.bbox).reshape(2,3).T for prop in regionprops_f])
    cent = np.array([np.around(prop.centroid).astype(int) for prop in regionprops_f])
    cell_id = np.arange(cent.shape[0])
    
    mask_plots = []
    mask_dims = []
    num_cells_list_im2 = []
    pcd_im2_list = []
    im2_patch_list = []    
    
    for i in (cell_id):
        bb = bboxes[i]
        bigmask = get_masks_with_edges(im2_bi>=1, bb, add_edge=add_edge)
        mask_plots.append(sparse.csr_matrix(bigmask.ravel()))
        mask_dims.append(bigmask.shape)        
    im2_mask = mask_plots
    im2_dims = mask_dims
    im2_cells = np.array([f.label for f in regionprops_f])
    im2_bb = bboxes
    
    for i in range(len(im2_mask)):
        img_i = np.array(im2_mask[i].todense())[0].reshape(im2_dims[i])
        numcells = cc3d.connected_components(img_i).max()
        num_cells_list_im2.append(numcells)
    
    for j in (range(len(im2_mask))):    
        img_j = np.array(im2_mask[j].todense())[0].reshape(im2_dims[j])
        pcd_j = o3d.geometry.PointCloud()
        pcd_j.points = o3d.utility.Vector3dVector(np.array(np.where(img_j)).T)                
        pcd_j_ds = pcd_j.voxel_down_sample(voxel_size=ds_sz)
        pcd_im2_list.append(pcd_j_ds)
    
    for j in (range(len(im2_mask))):            
        img = get_orig_2(orig_full=im2, bbox=im2_bb[j], add_edge=add_edge)
        if (np.array(img.shape)>=np.array(outputshape)).all():
            img_p = crop_w_given_sz(img, outputshape)
            img_p_z = zscore(img_p, axis = None)>2                
            img_p_ds = downscale_local_mean(img_p_z, (ds_sz_aff,ds_sz_aff,ds_sz_aff))        
            im2_patch_list.append(img_p_ds)    
        else:
            im2_patch_list.append(np.nan)    
    out = [im2_mask,im2_dims,im2_cells,im2_bb,num_cells_list_im2,pcd_im2_list,im2_patch_list]
    return(out)




def wahba(X,Y):
    X0=X-np.mean(X,axis=0)    
    Y0=Y-np.mean(Y,axis=0)    
    U, _, Vt = np.linalg.svd(X0.T@Y0)
    V = Vt.T
    M = np.eye(3)
    M[-1,:] = np.array((0,0,np.linalg.det(U)*np.linalg.det(V)))
    R = U@M@V.T
    T = np.mean(Y-X@R, axis = 0)
    return(np.r_[R,T[None,:]])  # return 4x3 matrix (to transform from first input to second inpupt)

def get_masks_with_edges(seg_cc_f, bbox, add_edge= np.array((10,100,100))):
    '''
    get the target position (specified by bbox) from seg_cc_f, 
    add edge too. 
    '''
    bb_with_edge = np.clip(bbox+np.c_[(-add_edge,add_edge)], np.zeros_like(bbox), np.c_[seg_cc_f.shape,seg_cc_f.shape]).astype(int)
    cen = bbox[:,0] + (np.ptp(bbox, axis = 1)/2)
    new_radius = np.min(np.abs(bb_with_edge-cen[:,None]), axis = 1)
    bb_with_edge_fixed = (np.c_[cen-new_radius, cen+new_radius]).astype(int)

    sls = [slice(x,y) for (x,y) in bb_with_edge_fixed]
    seg_cc_f = seg_cc_f[sls[0],sls[1],sls[2]]
    return(seg_cc_f)


def get_cmap_ordered(n, name = 'gist_ncar', grayscale = 0):
    '''
    n is the total number of neurons
    grayscale is the background, the higher the more white. 
    '''    
    import matplotlib
    np.random.seed(226)
    cmap = plt.cm.get_cmap(name, n)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    return(cmap)



def run_shape_icp(pcd_im1_list,num_cells_list_im1,
                  pcd_im2_list,num_cells_list_im2, 
                 thre = 3):
    '''
    im1 is the source, im2 is the target. 
    '''
    max_iteration = 100
    with_scaling=True
    transform_mtx_list = np.zeros((len(pcd_im2_list),len(pcd_im1_list), 4, 4))    
    all_eul = np.zeros((transform_mtx_list.shape[:2] + (3,))) + np.inf
    all_scl = np.zeros((transform_mtx_list.shape[:2])) + np.inf
    for i in (range(len(pcd_im2_list))):
        pcd_i_ds = pcd_im2_list[i]    
        for j in range(len(pcd_im1_list)):
            pcd_j_ds = pcd_im1_list[j]        
            if min(num_cells_list_im2[i], num_cells_list_im1[j])< 6:
                transform_mtx_list[i,j] = np.zeros((4,4))
                continue;        
            reg_p2p = o3d.pipelines.registration.registration_icp(pcd_j_ds,pcd_i_ds, max_correspondence_distance = thre, 
                                                                  estimation_method = o3d.pipelines.registration.TransformationEstimationPointToPoint(with_scaling=with_scaling), 
                                                                  criteria = o3d.cpu.pybind.pipelines.registration.ICPConvergenceCriteria(max_iteration=max_iteration),)
            scale = np.cbrt(np.linalg.det(reg_p2p.transformation[:3,:3]))
            if scale < 2 and scale > 0.6:
                transform_mtx_list[i,j] = reg_p2p.transformation            
                all_scl[i,j]=scale
                all_eul[i,j]=rotationMatrixToEulerAngles(reg_p2p.transformation[:3,:3]/scale)
            else:
                continue
    return(transform_mtx_list)
                
def run_tranform_get_NCC(transform_mtx_list,ds_sz_aff,
                         im1_patch_list,im2_patch_list):  
    '''
    again im1 is the source, im2 is the target. 
    '''
    icp_NCC = np.zeros(transform_mtx_list.shape[:2])    
    for i in (range(len(im2_patch_list))):
        im2p_ds = im2_patch_list[i]
        if np.isnan(im2p_ds).all():
            continue;    
        for j in range(len(im1_patch_list)):
            im1p_ds = im1_patch_list[j]
            if np.isnan(im1p_ds).all():
                continue;              
            else:            
                transform_mtx = np.array(transform_mtx_list[i,j])
                if (transform_mtx==0).all():
                    continue;
                offset = -(np.linalg.inv(transform_mtx)[:3,:3])@(transform_mtx[:-1, -1])/ds_sz_aff
                tranformed = ndi.affine_transform(im1p_ds, np.linalg.inv(transform_mtx)[:3,:3], offset=offset)          
                assert(im1p_ds.shape == im2p_ds.shape)  
                icp_NCC[i,j]=ncc(im2p_ds,tranformed)                            
    return(icp_NCC)



def find_consensus_parameters(icp_NCC,im1_cent,im2_cent,
                              im1_bi,im1_bb,im1_cells,
                              im2_bi,im2_bb,im2_cells,
                              add_edge=np.array((10,100,100)), T=50, max_num_cellpairs = 100):
    '''
    im1 is source, im2 is target
    '''
    topmatch_ids = np.array(np.unravel_index(np.argsort(icp_NCC.ravel())[::-1][:T], (icp_NCC.shape))).T    
#     ncc_thre = icp_NCC[topmatch_ids[-1,0], topmatch_ids[-1,1]]
    d_thre_list = np.arange(5,40,step=5)
    appearance_matrix_dict = {}
    for d_thre in d_thre_list:
        appearance_matrix_dict[d_thre] = np.zeros_like(icp_NCC, dtype = 'uint16')
    for k,(i,j)in (enumerate(topmatch_ids)):    
        foo1  = get_masks_with_edges(im2_bi,bbox=im2_bb[i], add_edge=add_edge)
        foo2  = get_masks_with_edges(im1_bi,bbox=im1_bb[j], add_edge=add_edge)
        foo1,foo2 = crop(foo1, foo2)    
        i_s = np.unique(foo1[foo1>0])  ## corresponds to the im2
        j_s = np.unique(foo2[foo2>0])  ## corresponds to the im1

        m_cent = np.array([(np.min(np.array(np.where(foo1 == m)).T, 0)+np.max(np.array(np.where(foo1 == m)).T, 0))/2 for m in i_s])
        n_cent = np.array([(np.min(np.array(np.where(foo2 == n)).T, 0)+np.max(np.array(np.where(foo2 == n)).T, 0))/2 for n in j_s])
        mn_d = distance.cdist(m_cent, n_cent)

        for d_thre in d_thre_list:
            for m,n in np.array(np.where(mn_d<d_thre)).T:
                i_p = np.where(im2_cells==i_s[m])[0][0]
                j_p = np.where(im1_cells==j_s[n])[0][0]
                appearance_matrix_dict[d_thre][i_p,j_p]+=1
    assert(np.array([mtx.max() for mtx in appearance_matrix_dict.values()]).max()>0)    
    rez = dict()
    for i,(d_thre,appearance_matrix) in enumerate(appearance_matrix_dict.items()):
        if appearance_matrix.max()==0:
            continue        
        a,b = np.unique(appearance_matrix, return_counts = True)        
        thre=1
        num_cells = 0
        for i in range(len(a)):
            if b[::-1][:i+1].sum() < max_num_cellpairs:
                thre = a[::-1][i] 
                num_cells = b[::-1][:i+1].sum()            
        pair_d = distance.cdist(im2_cent[np.where(appearance_matrix>=thre)[0]], im1_cent[np.where(appearance_matrix>=thre)[1]])
        cells_max_dist = np.diag(pair_d).max()
        rez[d_thre] = dict(thre = thre,num_cells = num_cells,cells_max_dist = cells_max_dist)
#     pprint(rez)
    final_d_thre = 0
    maxdist=np.inf
    final_thresh = 0
    for d_thre, v in rez.items():
        if (v['cells_max_dist']<1000) and v['num_cells']>=4:        
            if v['cells_max_dist']<maxdist:
                maxdist = v['cells_max_dist']
                final_d_thre = d_thre
                final_thresh = v['thre']
                appearance_matrix = appearance_matrix_dict[d_thre]        
        else:
            continue;
    print(f'final distance threshold = {final_d_thre}\nfinal coappearance threshold = {final_thresh}')        
    
    
    appearance_matrix_orig = appearance_matrix.copy()
    appearance_matrix = rm_one2many_from_app_mtx(appearance_matrix, final_thresh)

    return(final_d_thre, final_thresh, appearance_matrix)



def get_cellcent(im2_bi):
    cols = ['z','x','y']
    exvivo_centroids = np.array([region.centroid for region in measure.regionprops(im2_bi)])
    df_im2 = (pd.DataFrame(exvivo_centroids, columns = ['z','x','y']))[['x','y','z']]
    df_im2 = df_im2[cols]
    im2_cent = df_im2.values
    return(im2_cent)


def save_MRS_file(max_rot_set_f, im1_cent,im2_cent,appearance_matrix,final_thresh):
    X = im2_cent[np.array(np.where(appearance_matrix>=final_thresh))[0]]
    Y = im1_cent[np.array(np.where(appearance_matrix>=final_thresh))[1]]
    if X.shape[0]>20:
        X =X[:20]
        Y =Y[:20]
    X_dst = distance.squareform(distance.pdist(X))
    Y_dst = distance.squareform(distance.pdist(Y))
    # max_rot_set_f = f'/home/ubuntu/immunostaining/simulation/MRS/'
    if not os.path.exists(max_rot_set_f):
        os.mvkdir(max_rot_set_f)    
    savemat(max_rot_set_f + f'dst_X.mat', {'X_dst':X_dst})
    savemat(max_rot_set_f + f'dst_Y.mat', {'Y_dst':Y_dst})    
    
def get_exact_matching_from_MRS(max_rot_set_f, im1_cent, im2_cent, appearance_matrix,final_thresh):
    '''
    after MRS, return bhat. 
    '''
    rot_dist_lb = loadmat(max_rot_set_f + f'rot_dist_lb.mat')
    rot_dist_lb = rot_dist_lb['rot_dist_lb']
    final_set = np.where(np.diag(rot_dist_lb)==rot_dist_lb.max())[0]
    final_set_id = np.array(np.where(appearance_matrix>=final_thresh)).T[final_set]
    X_p=im2_cent[final_set_id[:,0]]
    Y_p=im1_cent[final_set_id[:,1]]
    bhat = wahba(Y_p,X_p)  # transformm from Y to X
    return(bhat, X_p,Y_p)



def get_neighbor_coord(im1_cent,im2_cent, im2, X_p,Y_p):
    brdr=np.array((20,200,200))
#     brdr=80
    mn_im2 = np.clip(np.min(X_p, 0).astype(int)-brdr, 0, None)
    mx_im2 = np.clip(np.max(X_p, 0).astype(int)+brdr, None, im2.shape)
    im2_consensuscell_ids = np.where([((l>=mn_im2).all() and (l<=mx_im2).all()) for l in im2_cent])[0]

    mn_im1=np.clip(np.min(Y_p, 0) - brdr, 0, None)
    mx_im1=np.max(Y_p, 0) + brdr
    im1_consensuscell_ids = np.where([((l>=mn_im1).all() and (l<=mx_im1).all()) for l in im1_cent])[0]

    set1p_cent = im1_cent[im1_consensuscell_ids]
    set2p_cent = im2_cent[im2_consensuscell_ids]
    return(set1p_cent, set2p_cent)


def run_doubleICP(im1, im2, bhat, set1p_cent, set2p_cent):
    tran_shape=np.array(im2.shape)
    max_correspondence_distance = 20
    max_iteration = 100
    foo = np.c_[set1p_cent, np.ones((set1p_cent.shape[0],1))]   #Nx4
    set1p_cent_t = foo@bhat
    
    offset = -(bhat[:3,:3])@bhat[-1]
    im1_wahba = ndi.affine_transform(im1, bhat[:3,:3],offset = offset, output_shape = tran_shape)    

    pcd1 = o3d.geometry.PointCloud()
    pcd1.points = o3d.utility.Vector3dVector(set1p_cent_t)
    pcd2 = o3d.geometry.PointCloud()
    pcd2.points = o3d.utility.Vector3dVector(set2p_cent)
    reg_p2p_1 = o3d.pipelines.registration.registration_icp(pcd1,pcd2, max_correspondence_distance = max_correspondence_distance, 
                                                          estimation_method = o3d.pipelines.registration.TransformationEstimationPointToPoint(with_scaling=True), 
                                                          criteria = o3d.cpu.pybind.pipelines.registration.ICPConvergenceCriteria(max_iteration=max_iteration),)
    tran = reg_p2p_1.transformation
    
    while np.isnan(tran).sum()>0:
        max_correspondence_distance += 10
        reg_p2p_1 = o3d.pipelines.registration.registration_icp(pcd1,pcd2, max_correspondence_distance = max_correspondence_distance, 
                                                          estimation_method = o3d.pipelines.registration.TransformationEstimationPointToPoint(with_scaling=True), 
                                                          criteria = o3d.cpu.pybind.pipelines.registration.ICPConvergenceCriteria(max_iteration=max_iteration),)
        tran = reg_p2p_1.transformation
    return(tran, im1_wahba)


def get_final_transformation(tran,bhat):
    A = tran[:3,:3]
    t_A = tran[:3,-1]
    B = bhat[:3,:3]
    t_B = bhat[-1]
    total_shift = t_B@A.T + t_A
    total_R = (B@A.T)  # i dont know why theres an T in the end but this works! - removed. 
    total_scale = np.cbrt(np.linalg.det(total_R))
    out = [total_R, total_shift, total_scale]
    return(out)


def apply_total_R(im1,im2,total_R, total_shift,):
    t = -total_shift.T@np.linalg.inv(total_R)
    im1_happy = ndi.affine_transform(im1, np.linalg.inv(total_R.T), offset = t, output_shape=im2.shape)    
    return(im1_happy)


def apply_doubleICP_2_binary(im1_bi,im2_bi,total_R,total_shift):
    im2_bi_01 = (im2_bi>0).astype(bool)
    im1_bi_01 = (im1_bi>0).astype(bool)
    t = -total_shift.T@np.linalg.inv(total_R)
    im1_bi_happy = ndi.affine_transform(im1_bi_01, np.linalg.inv(total_R.T), offset = t, output_shape=im2_bi.shape, order = 0)        
    return(im1_bi_happy)
    
def apply_deformable_2_binary(im1_happy,im1_bi_happy,im1_bi,im2_bi,vec_ds,vec_field_smooth_list):    
    im2_bi_01 = im2_bi>0
    im1_bi_01 = im1_bi>0    
    im2_bi_dsz = downscale_local_mean(im2_bi_01, (vec_ds,vec_ds,vec_ds))
    im1_bi_doubleicp_dsz = downscale_local_mean(im1_bi_happy, (vec_ds,vec_ds,vec_ds)).astype(bool)
    out = im1_bi_doubleicp_dsz.copy()  # out is ds-ed
    ncc_0 = np.around(ncc(im2_bi_dsz,im1_bi_doubleicp_dsz),3)
    for i in range(len(vec_field_smooth_list)): 
#         if i%5==0:
#             print(f'iteration {i}...')
        vec_field_smooth = vec_field_smooth_list[i]
        assert((im1_bi_doubleicp_dsz.shape == vec_field_smooth.shape[:-1]))
        shape = im1_bi_doubleicp_dsz.shape
        mapz_base, mapx_base, mapy_base = np.meshgrid(np.arange(shape[0]),np.arange(shape[1]),np.arange(shape[2]), indexing='ij')
        mapz=mapz_base-vec_field_smooth[:,:,:,0]
        mapx=mapx_base-vec_field_smooth[:,:,:,1]
        mapy=mapy_base-vec_field_smooth[:,:,:,2]
        im1_bi_de = warp(out,np.array((mapz,mapx,mapy)), order = 0)
        out = im1_bi_de
    im1_bi_de = rescale(out,(vec_ds,vec_ds,vec_ds), order = 0)
    return(im1_bi_de)



def count_cell_matching(im1_bi_de, im2_bi):
    im2_bi_01 = im2_bi>0
    im1_bi_de_01 = im1_bi_de.astype(bool)
    cols = ['z','x','y']
    cc_im2 = cc3d.connected_components(im2_bi_01)
    foo = np.array([region.centroid for region in measure.regionprops(cc_im2)])
    df_im2_new = (pd.DataFrame(foo, columns = ['z','x','y']))[['x','y','z']]
    df_im2_new = df_im2_new[cols]
    im2_cent_new = df_im2_new.values

    cc_im1 = cc3d.connected_components(im1_bi_de_01)
    foo = np.array([region.centroid for region in measure.regionprops(cc_im1)])
    df_im1_new = (pd.DataFrame(foo, columns = ['z','x','y']))[['x','y','z']]
    df_im1_new = df_im1_new[cols]
    im1_cent_new = df_im1_new.values
    dst = distance.cdist(im2_cent_new, im1_cent_new)    
    closest_pair =[]
    for j,i in enumerate(np.argmin(dst,0)):
        vec = dst[i]
        if np.argmin(vec)==j:
            if dst[i,j] <= 30:
                closest_pair.append((i,j))
    closest_pair = np.array(closest_pair)    
    print('***FOUND PAIRS:***',closest_pair.shape[0])
    out = [cc_im1,cc_im2,closest_pair,im1_cent_new]
    return(out)


def rm_one2many_from_app_mtx(appearance_matrix, final_thresh):    
    out=appearance_matrix.copy()
    cellpair = np.array(np.where(appearance_matrix>=final_thresh)).T    
    a = cellpair[:,0]
    unq, unq_idx, unq_cnt = np.unique(a, return_inverse=True, return_counts=True)
    cnt_mask = unq_cnt > 1
    dup_entry = unq[cnt_mask]
    dup_ids = np.where([(m in dup_entry) for m in a])[0]
    for p in cellpair[dup_ids]:
        out[p[0],p[1]]=0
    cellpair = np.delete(cellpair, dup_ids, axis=0)

    a = cellpair[:,1]
    unq, unq_idx, unq_cnt = np.unique(a, return_inverse=True, return_counts=True)
    cnt_mask = unq_cnt > 1
    dup_ids = unq[cnt_mask]
    dup_entry = unq[cnt_mask]
    dup_ids = np.where([(m in dup_entry) for m in a])[0]
    for p in cellpair[dup_ids]:
        out[p[0],p[1]]=0
    cellpair = np.delete(cellpair, dup_ids, axis=0)
    return(out)    





def gray2RGB(img, col, scl = 1):
    cols = ['r','g','b','y','m','cyan','white','purple','pink', 'orange']
    if len(img.shape)==3:
        img = np.max(img,axis=0)
    assert(len(img.shape)==2)
    colid = np.where(np.array([col==c for c in cols]))[0][0]    
    img_rgb = np.zeros(img.shape+(3,))
    if colid <3: # RGB case:        
        img_rgb[:,:,colid] = img.copy()
    else:
        if col=='y': # if its yellow, 
            img_rgb[:,:,0] = img.copy()
            img_rgb[:,:,1] = img.copy()
            img_rgb[:,:,2] = img.copy()*0
        if col=='m': # magenta case,
            img_rgb[:,:,0] = img.copy()
            img_rgb[:,:,2] = img.copy()
        if col=='cyan': # cyan
            img_rgb[:,:,1] = img.copy()
            img_rgb[:,:,2] = img.copy()
        if col=='white': # white
            img_rgb[:,:,0] = img.copy()
            img_rgb[:,:,1] = img.copy()
            img_rgb[:,:,2] = img.copy()
        if col=='purple': # purple
            img_rgb[:,:,0] = img.copy()*0.33
            img_rgb[:,:,1] = img.copy()*0.33
            img_rgb[:,:,2] = img.copy()*0.80
        if col=='pink': # pink
            img_rgb[:,:,0] = img.copy()
            img_rgb[:,:,1] = img.copy()*0.41
            img_rgb[:,:,2] = img.copy()*0.70
        if col=='orange': # pink
            img_rgb[:,:,0] = img.copy()
            img_rgb[:,:,1] = img.copy()*0.784
            img_rgb[:,:,2] = img.copy()*0
    img_rgb /=img_rgb.max()*scl
    return(img_rgb)

