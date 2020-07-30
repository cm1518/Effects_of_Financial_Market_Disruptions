function [b0 stdl0 res0 sigma]=VAR_SR(Y)
global l 

%we drop the first l infos bec of the lag
T=size(Y,1);
dim=size(Y,2);

Y0=Y(l+1:T,:);
t=(l:T-1);

for i=1:l
Y_{i}=Y((l+1-i):(T-i),:);
end

%WE NOW COMPUTE THE UNRESTRICTED VAR (NO QUARTERLY DEPENDENCE)
X0=[];
for i=1:l
    X0=[X0 Y_{i}];
end;
X0=[X0 ones(T-l,1)];

b0=inv(X0'*X0)*X0'*Y0;
A_OLS = inv(X0'*X0)*(X0'*Y0);
SSE = (Y0 - X0*A_OLS)'*(Y0 - X0*A_OLS);
sigma = SSE./(length(Y0)-length(b0));
% [bb_y,bint,r,rint,stats]=  regress(y,X0)
% eg: b0=[b0_y b0_pi b0_r];

res0=Y0-X0*b0;
stdl0=0;

