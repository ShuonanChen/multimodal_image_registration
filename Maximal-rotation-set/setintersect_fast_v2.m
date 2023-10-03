function [out1,out2]=keepintersect_fast_v2(a,b)
% Credit - Divakar @ stackoverflow
%// Calculate equivalent one-column versions of input arrays
mult = (10^ceil(log10( 1+max( [a(:);b(:)] ))).^(size(a,2)-1:-1:0))'; %//'
acol1 = a*mult;
bcol1 = b*mult;

%// Use ismember to get indices of the common elements
[match_a,~] = ismember(acol1,bcol1);

%// Now, with ismember, duplicate items are not taken care of automatically as
%// are done with intersect. So, we need to find the duplicate items and
%// remove those from the outputs of ismember
[~,a_sorted_ind] = sort(acol1);
a_rm_ind =a_sorted_ind([false;diff(sort(acol1))==0]); %//indices to be removed
match_a(a_rm_ind)=0;

ind_a = find(match_a);

out1 = a(ind_a,1);
out2 = a(ind_a,2);

end