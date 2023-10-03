function [out1,out2]=mykeepintersect(x,y)
a=ismember(x,y,'rows');
out1=x(a==1,1);
out2=x(a==1,2);
end