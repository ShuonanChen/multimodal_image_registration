function [R,T,beta]=wahba(X,Y)
% Wahba's algorithm for solving for the least-squares rotation
% returns solution to ||Y - X*R - T||, R \in S0(3)
% alternate form ||Y - [X 1]*beta||


X0=X-mean(X,1);
Y0=Y-mean(Y,1);

[U,~,V]=svd(X0'*Y0);

% M=[1 0 0;0 1 0;0 0 det(U)*det(V)];
M=eye(size(X,2));
M(end,end)=det(U)*det(V);
R=U*M*V';


T= mean(Y - X*R,1);

beta = [R;T];
end