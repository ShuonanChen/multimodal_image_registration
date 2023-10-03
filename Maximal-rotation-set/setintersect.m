function [out1,out2]=setintersect(x,y)
% code to find the common rows between two matrices. There are 3 routines,
% one is considerably faster than others thanks to fast matrix mult. algs.

% [out1,out2]=setintersect_fast_v1(x,y);
[out1,out2]=setintersect_fast_v2(x,y);
% [out1,out2]=mysetintersect(x,y);
end