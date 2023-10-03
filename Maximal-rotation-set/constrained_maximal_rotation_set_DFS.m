function [rot_set_lower_bound,rot_set_upper_bound]=constrained_maximal_rotation_set_DFS(DX,DY,margin,maxdepth,X,Y,anglelimit,include)
%Globally optimal algorithm to find the maximal number of pairs of points
%in two datasets that are rotations of one another (up to a margin of
%error). ** Rotations are constrained to reside in a particular angular
%range.
%Inputs - DX,DY - are nxn and mxm pairwise distance matrices
% margin - the error tolerance for the matching |DX(i,j) - DY(k,l)|<margin
% to be considered a match
% maxdepth - maximum depth for Branch and Bound -- if maxdepth>=min(m,n),
% this algorithm will return the global optima
% X,Y - are n x 3 and m x 3 spatial locations of points that are used to
% compute angles.
% anglelimit - 1 x 3 vector that constrains the maximum euler angles for
% the rotation (in degrees)
% include - binary n x m matrix to constrain which pairs (i,k) between i in X and k in Y are
% included in the matching (if empty, include everything)
%Outputs - rot_set_lower_bound - lower bound on the maximal set of pairs of
%points in X and Y such that X_i and Y_j are in it
% rot_set_upper_bound - upper bound on the maximal set of pairs of
%points in X and Y such that X_i and Y_j are in it
%
% Note that if lower bound and upper bound are equal, this is the maximal
% set size for pairs of points X_i and Y_j.
% If all lower bounds equal all upper bounds, then maximal lower bound is
% the global optimum maximal set size.
if nargin<9
    include=ones(size(DX,1),size(DY,1));
end
x_idx=cell(size(DX,1),size(DY,1));
y_idx=cell(size(DX,1),size(DY,1));
global rot_set_lower_bound rot_set_upper_bound%#ok<REDEFGI>
%initial graph building
for i=1:size(DX,1)
    for j=1:size(DY,1)
        [II0,JJ0]=find(abs(DX(i,:)'-DY(j,:))<margin);
        keep=find((II0~=i).*(JJ0~=j));
        x_idx{i,j}=II0(keep);
        y_idx{i,j}=JJ0(keep);
    end
end
rot_set_lower_bound=ones(size(DX,1),size(DY,1));
rot_set_upper_bound=ones(size(DX,1),size(DY,1));
tic;t=0;numsearch=size(DX,1)*size(DY,1);
for i=1:size(DX,1)
    for j=1:size(DY,1)
        
        t=t+1;
        if include(i,j)==1
            II=x_idx{i,j};JJ=y_idx{i,j};
            recursion(II,JJ,x_idx,y_idx,1,maxdepth,i,j,[i],[j],DX,DY,margin,X,Y,anglelimit); %#ok<NBRAK> %recurse with (i,j)
        end
        clc
        fprintf(['Searching for transformations (' num2str(t) '/' num2str(numsearch) ')...\n']);
        fprintf(['\n' repmat('.',1,50) '\n\n'])
        for tt=1:round(t*50/(numsearch))
            fprintf('\b|\n');
        end
        TT=toc;
        disp(['Time elapsed (minutes): ' num2str(TT/60) ' Time remaining (minutes): ' num2str((numsearch-t)*(TT/t)*(1/60)) ' Est. Total (minutes): ' num2str(TT/60 + (numsearch-t)*(TT/t)*(1/60))]);
    end
end
end

% Check if a pair of point sets is a rotation
function out=isrotation(DX,DY,I,J,margin)
vec=@(x)(x(:));
if all(vec(abs(DX(I,I)-DY(J,J)))<margin)
    out=1;
else
    out=0;
end
end

% Check euler angles of a rotation (only if there are at least 4 points)
function out=angle_check(X,Y,cur_set_x,cur_set_y,anglelimit)
if length(cur_set_x)>3
    angles=rad2deg(rotm2eul(wahba(X(cur_set_x,:),Y(cur_set_y,:))));
    if all(abs(angles)<=anglelimit)
        out=1;
    else
        out=0;
    end
else
    out=1;
end

end

% Recursive function to branch and bound
function recursion(II,JJ,x_idx,y_idx,depth,maxdepth,i,j,cur_set_x,cur_set_y,DX,DY,margin,X,Y,anglelimit)
if angle_check(X,Y,cur_set_x,cur_set_y,anglelimit)==0 %% might need to put this lower at some place - has slight bugs
    return 
end

global rot_set_lower_bound rot_set_upper_bound

if or(or(depth==maxdepth,isempty(II)),length(II)<=rot_set_lower_bound(i,j)-length(cur_set_x)) % terminate if reached max depth or there is nothing to recurse on or if best future size is lower than current lower bound
    rot_set_upper_bound(i,j)=max([rot_set_upper_bound(i,j) length(cur_set_x)+length(II) rot_set_lower_bound(i,j)]); % switching this and lower bound orders might solve the bug
    rot_set_lower_bound(i,j)=max(rot_set_lower_bound(i,j),length(cur_set_x));
    return
end

if isrotation(DX,DY,cur_set_x,cur_set_y,margin)==1
    ind=sub2ind(size(rot_set_lower_bound),cur_set_x,cur_set_y);
    rot_set_lower_bound(ind)=max(rot_set_lower_bound(ind),length(cur_set_x));
end

if and(isrotation(DX,DY,[cur_set_x';II],[cur_set_y';JJ],margin)==1,angle_check(X,Y,[cur_set_x';II],[cur_set_y';JJ],anglelimit)==1)
    ind=sub2ind(size(rot_set_lower_bound),[cur_set_x';II],[cur_set_y';JJ]);
    rot_set_lower_bound(ind)=max(rot_set_lower_bound(ind),length(cur_set_x)+length(II));
end

for u=1:size(II,1)
    ii=II(u);jj=JJ(u);
    [III,JJJ]=setintersect([II JJ],[x_idx{ii,jj} y_idx{ii,jj}]);
    recursion(III,JJJ,x_idx,y_idx,depth+1,maxdepth,i,j,[cur_set_x ii],[cur_set_y jj],DX,DY,margin,X,Y,anglelimit);
end
end