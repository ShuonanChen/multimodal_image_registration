clear
clc
close all
cd('/home/ubuntu/code/Maximal-rotation-set')

load('./dst_X.mat')
load('./dst_Y.mat')
size(X_dst)

% creat mask
include = eye(length(X_dst))

resolution=20;
maxdepth=10;
tic;[rot_dist_lb,rot_dist_ub]=maximal_rotation_set_DFS(X_dst,Y_dst,resolution,maxdepth,include);DFS_toc=toc;
imagesc([rot_dist_lb]);
axis equal;axis tight
colorbar

maximum = max(max(rot_dist_lb));
[x,y]=find(rot_dist_lb==maximum)

% now x and y are the row and col indices of the max positions!
i0=x(1)
j0=y(1)
[bestset]=maximal_rotation_set_DFS_decode(X_dst,Y_dst,resolution,maxdepth,i0,j0)



save('./rot_dist_lb.mat','rot_dist_lb')

