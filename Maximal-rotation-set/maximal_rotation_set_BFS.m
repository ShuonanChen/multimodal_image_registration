function [rot_set_lower_bound,rot_set_upper_bound]=maximal_rotation_set_BFS(DX,DY,margin)
%Globally optimal algorithm to find the maximal number of pairs of points
%in two datasets that are rotations of one another (up to a margin of
%error). This algorithm performs breadth first search using DFS
%subroutines.
%Inputs - DX,DY - are nxn and mxm pairwise distance matriceces
% margin - the error tolerance for the matching |DX(i,j) - DY(k,l)|<margin
% to be considered a match
%Outputs - rot_set_lower_bound - lower bound on the maximal set of pairs of
%points in X and Y such that X_i and Y_j are in it
% rot_set_upper_bound - upper bound on the maximal set of pairs of
%points in X and Y such that X_i and Y_j are in it
%
% Note that if lower bound and upper bound are equal, this is the maximal
% set size for pairs of points X_i and Y_j.
% If all lower bounds equal all upper bounds, then maximal lower bound is
% the global optimum maximal set size.
include=ones(size(DX,1),size(DY,1));
finished=0;
maxdepth=0;
rot_set_lower_bound=zeros(size(include));
rot_set_upper_bound=inf(size(include));
while finished==0
    maxdepth=maxdepth+1;
    [dfs_lb,dfs_ub]=maximal_rotation_set_DFS(DX,DY,margin,maxdepth,include);
    rot_set_lower_bound(include==1)=max(dfs_lb(include==1),rot_set_lower_bound(include==1));
    rot_set_upper_bound(include==1)=min(dfs_ub(include==1),rot_set_upper_bound(include==1));
    include=(dfs_ub>dfs_lb);
    if all(include(:)==0)
        finished=1;
    end
end



end