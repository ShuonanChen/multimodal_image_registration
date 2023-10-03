function [bestset]=maximal_rotation_set_DFS_decode(DX,DY,margin,maxdepth,i0,j0);

x_idx=cell(size(DX,1),size(DY,1));
y_idx=cell(size(DX,1),size(DY,1));
global rot_set_lower_bound rot_set_upper_bound bestset %#ok<REDEFGI>
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
bestset=[];

II=x_idx{i0,j0};JJ=y_idx{i0,j0};
recursion(II,JJ,x_idx,y_idx,1,maxdepth,i0,j0,[i0],[j0],DX,DY,margin); %#ok<NBRAK> %recurse with (i,j)
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

% Recursive function to branch and bound
function recursion(II,JJ,x_idx,y_idx,depth,maxdepth,i,j,cur_set_x,cur_set_y,DX,DY,margin)
global rot_set_lower_bound rot_set_upper_bound bestset

if or(or(depth==maxdepth,isempty(II)),length(II)<=rot_set_lower_bound(i,j)-length(cur_set_x)) % terminate if reached max depth or there is nothing to recurse on or if best future size is lower than current lower bound
    rot_set_upper_bound(i,j)=max([rot_set_upper_bound(i,j) length(cur_set_x)+length(II) rot_set_lower_bound(i,j)]);
    rot_set_lower_bound(i,j)=max(rot_set_lower_bound(i,j),length(cur_set_x));
    return
end

if isrotation(DX,DY,[cur_set_x';II],[cur_set_y';JJ],margin)==1
    ind=sub2ind(size(rot_set_lower_bound),[cur_set_x';II],[cur_set_y';JJ]);
    if rot_set_lower_bound(i,j)<length(cur_set_x)+length(II)
        bestset=[[cur_set_x';II] [cur_set_y';JJ]];
    end
    rot_set_lower_bound(ind)=max(rot_set_lower_bound(ind),length(cur_set_x)+length(II));
end

for u=1:size(II,1)
    ii=II(u);jj=JJ(u);
    [III,JJJ]=setintersect([II JJ],[x_idx{ii,jj} y_idx{ii,jj}]);
    recursion(III,JJJ,x_idx,y_idx,depth+1,maxdepth,i,j,[cur_set_x ii],[cur_set_y jj],DX,DY,margin);
end
end