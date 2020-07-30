close all

global gdpm pgdpm pcom fyff iter A Ap An A0 Om_1 l n var Q1 Q2 Q3 Q4 psi0 PSI Om0 Om1 Om2 Om3 Om4 s trend W U0 indice

%XXXXXXXXXXXXXXXXXXXXXXXXXXXX     INPUT    XXXXXXXXXXXXXXXXXXXXXXXXXXXX

% l=number of lags
l=setup.VARsym_order;
%number of lags for IRF
n=setup.lags;
%number of lags for MA(infinie)
s=setup.lags;
%size of bootstrap
nt=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        EMPIRICAL SPECIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFICATION : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load gzdata
gdpm = data2(:,1); pgdpm = data2(:,2); pcom = data2(:,3); fyff = data2(:,4); 
variablename={'GDP Growth','Inlation','Excess Bond Premium','Fed Funds Rate'};
Y=[gdpm pgdpm pcom fyff];
%get rid of missing (NaN) observations
[i j]=find((isnan(Y)==1));
Y(1:max(i),:)=[];

gdpm=Y(:,1);
pgdpm=Y(:,2);
pcom=Y(:,3);
fyff=Y(:,4);

data=Y'; %use for storage later
figure, plot(Y)
% varname={'Def','NonDef','Tax','GDP'};
legend(variablename)
% I find the VAR correponding to the symmetric case
[b0 stdl0 res0 sigma]= VAR_SR(Y);
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
A=b0(1:end-1,:)'; %I remove the constant coefs to fit in CM's code
data=[gdpm pgdpm pcom fyff]';
Om=cov(res0);
Om_1=inv(Om);
Ai0=chol(Om)';

%% Fully Recursive FAIR
A0=Ai0;

%% In case of proxy FAIR
load A_psvar
Ai0 = A_full;
A0 = A_full;

var=size(b0,2);
C=MA0(b0,res0,s);

A=[];
for ju=2:size(C,3)
    A(:,:,ju-1)=C(:,:,ju)*A0;
end

eps=inv(A0)*res0';

N=length(eps);Y=[];
for j=1:N
    Y(:,j)=eps(:,j);
    for k=1:min(length(A),j-1)
        Ah=A(:,:,k);
        Y(:,j)=Y(:,j)+Ah*eps(:,j-k);
    end
end
Y=Y';

Ap=A;
An=A;


IRF_mp(:,1)=[A0(1,3); reshape(A(1,3,1:s),1,s)']; %Y
IRF_mp(:,2)=[A0(2,3); reshape(A(2,3,1:s),1,s)']; %Pi
IRF_mp(:,3)=[A0(3,3); reshape(A(3,3,1:s),1,s)']; %EBP
IRF_mp(:,4)=[A0(4,3); reshape(A(4,3,1:s),1,s)']; %FFr


figure, plot(IRF_mp)
legend(variablename(1:end))
title('Impulse Responses to EBP Shock')

save ((['data_file',num2str(setup.size_obs),'.mat']), 'data')


