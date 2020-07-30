function [psi0]=MA0(b0, res0,s)
%s=number of coef in MA

global var l Om0

% n: length of IRF
% b:matrix of estimated coeficients

phi=zeros(var,var,l);
psi0=zeros(var,var,s+1);

for j=1:l
    phi(:,:,j)=b0((j-1)*var+1:(j*var),:)';
end;

psi0(:,:,1)=eye(var,var);

for j=1:s
    for i=1:l
        if j+1-i>0
        psi0(:,:,1+j)=phi(:,:,i)*psi0(:,:,1+j-i)+psi0(:,:,1+j);
        end;
    end;
end;