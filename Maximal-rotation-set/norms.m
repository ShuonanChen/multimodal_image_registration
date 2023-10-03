function n=norms(x,dim)
% row-wise norms
n=sqrt(sum(x.^2,dim));
end
