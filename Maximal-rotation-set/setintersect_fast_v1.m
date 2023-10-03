function [out1,out2]=keepintersect_fast_v1(a,b)
% Credit - Divakar @ stackoverflow
%// Calculate equivalent one-column versions of input arrays
mult = (10^ceil(log10( 1+max( [a(:);b(:)] ))).^(size(a,2)-1:-1:0))'; %//'
acol1 = a*mult;
bcol1 = b*mult;

[~, ind_a, ~] = intersect(acol1,bcol1);
out1 = a(ind_a,1);
out2 = a(ind_a,2);

end