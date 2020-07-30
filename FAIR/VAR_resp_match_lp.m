%clear all
%clc

%load x_new
%setup_structure_flexible_opt;
setup.horizon=setup.lags; %horizon up to which IRFS are matched

for kk=1:setup.horizon+1
setup.store_responses(:,:,kk)=VAR_MA_implied_higher_order( setup,kk-1);
end


%making the code backwards compatible
%load irf3
%setup.store_responses(3,3,1:61) = irf3;
%setup.store_responses(3,3,62:end) = 0;
store_responses=setup.store_responses;



% 1.0515

options = optimset('Display','iter','TolFun',1e-20,'TolX',1e-20,'MaxIter',5000, 'MaxFunEvals',500000);



%first optimize for only 1 Gaussian

for j=1:setup.size_obs
    for k=1:setup.size_obs
    if k>j %above the diagonal
        opt_restr=@(params)var_quad_restr( params,setup,store_responses,j,k ) ;
%        [xestimate,functionvalue1]=fminsearch(opt_restr,ones(3,1),options);
if j == 1;
       [xestimate,functionvalue1]=fminsearch(opt_restr,1*[-0.014;3; 5],options);
else 
           [xestimate,functionvalue1]=fminsearch(opt_restr,1*[-0.3;20; 1000],options);
         %[xestimate,functionvalue1]=fmincon(opt_restr,[-0.3;9; 20],[],[],[],[],[-1000,0.1,-1000]',[1000,8,100]',[],options);

end
        contemp(j,k)=0;
        beta_gen(j,k)=xestimate(1);
        b_gen(j,k)=xestimate(2);
        c_gen(j,k)=xestimate(3);
elseif j==k
        opt_free=@(params)var_quad_free( params,setup,store_responses,j,k ) ;
%        [xestimate,functionvalue1]=fminsearch(opt_free,ones(4,1),options);
       [xestimate,functionvalue1]=fminsearch(opt_free,2*ones(4,1),options);
        % [xestimate, functionvalue1]= fmincon(opt_free,[1,1,1,1],[],[],[],[],[-1000; -1000; 0.01; -1000],[],[],options)

        
        beta_diag(j)=(xestimate(1));
        
        beta_gen(j,k)=xestimate(2);
        b_gen(j,k)=xestimate(3);
        c_gen(j,k)=xestimate(4);
    elseif k<j
    opt_below_diag=@(params)var_quad_below_diag( params,setup,store_responses,j,k ) ;
%        [xestimate,functionvalue1]=fminsearch(opt_free,ones(4,1),options);
      [xestimate,functionvalue1]=fminsearch(opt_below_diag,5*ones(3,1),options);
      %           [xestimate,functionvalue1]=fmincon(opt_below_diag,[-0.3;9; 20],[],[],[],[],[-1000,0.1,-1000]',[1000,25,2200]',[],options);

       
        
        beta_gen(j,k)=xestimate(1);
        b_gen(j,k)=xestimate(2);
        c_gen(j,k)=xestimate(3);
    end
    
    end
end


%then optimize for two Gaussians

for j=1:setup.size_obs
    for k=1:setup.size_obs
    if k>j %above the diagonal
        opt_restr=@(params)var_quad_restr_2_gaussian( params,setup,store_responses,j,k ) ;
%        [xestimate,functionvalue1]=fminsearch(opt_restr,ones(3,1),options);
if j == 1;
        %   [xestimate, functionvalue1]= fmincon(opt_restr,[beta_gen(j,k); b_gen(j,k);c_gen(j,k); -beta_gen(j,k); 5*b_gen(j,k);2*c_gen(j,k)],[],[],[],[],[-1000; 0.1; -1000; -1000; 11; -1000],[1000; 4; 6.5; 1000; 14; 5 ],[],options)
       %    [xestimate, functionvalue1]= fmincon(opt_restr,[beta_gen(j,k); b_gen(j,k);c_gen(j,k); -beta_gen(j,k); 2*b_gen(j,k); 5*c_gen(j,k)],[],[],[],[],[-1000; 0.1; -1000; -1000; 32; -1000],[],[],options)

       [xestimate,functionvalue1]=fminsearch(opt_restr,[beta_gen(j,k); b_gen(j,k);c_gen(j,k);0.1; 35;30],options);
else
        %   [xestimate,functionvalue1]=fminsearch(opt_restr,[beta_gen(j,k); b_gen(j,k);c_gen(j,k);beta_gen(j,k); 2*b_gen(j,k);c_gen(j,k)],options);
        [xestimate, functionvalue1]= fmincon(opt_restr,[beta_gen(j,k); b_gen(j,k);c_gen(j,k); beta_gen(j,k); b_gen(j,k);c_gen(j,k)],[],[],[],[],[-1000; 10; -1000; -1000; 30; -1000],[0; 20; 1000; 0; 40; 1000],[],options)


end
       contemp12(j,k)=0;
        beta_gen12(j,k)=xestimate(1);
        b_gen12(j,k)=xestimate(2);
        c_gen12(j,k)=xestimate(3);
        beta_gen2(j,k)=xestimate(4);
        b_gen2(j,k)=xestimate(5);
        c_gen2(j,k)=xestimate(6);
    elseif j==k
        opt_free=@(params)var_quad_free_2_gaussian( params,setup,store_responses,j,k ) ;
%        [xestimate,functionvalue1]=fminsearch(opt_free,ones(4,1),options);
        %[xestimate,functionvalue1]=fminsearch(opt_free,[beta_diag(j);beta_gen(j,k); b_gen(j,k);c_gen(j,k);beta_gen(j,k); 10*b_gen(j,k);c_gen(j,k)],options);
        [xestimate, functionvalue1]= fmincon(opt_free,[beta_diag(j);beta_gen(j,k); b_gen(j,k);c_gen(j,k);beta_gen(j,k); 10*b_gen(j,k);c_gen(j,k)],[],[],[],[],[-1000; -1000; 0.1; -1000; -1000; 10; -1000],[],[],options)
        beta_diag12(j)=(xestimate(1));
        
        beta_gen12(j,k)=xestimate(2);
        b_gen12(j,k)=xestimate(3);
        c_gen12(j,k)=xestimate(4);
        beta_gen2(j,k)=xestimate(5);
        b_gen2(j,k)=xestimate(6);
        c_gen2(j,k)=xestimate(7);
    elseif k<j
    opt_below_diag=@(params)var_quad_below_diag_2_gaussian( params,setup,store_responses,j,k ) ;
%        [xestimate,functionvalue1]=fminsearch(opt_free,ones(4,1),options);
        [xestimate,functionvalue1]=fminsearch(opt_below_diag,[beta_gen(j,k); b_gen(j,k);c_gen(j,k);beta_gen(j,k); b_gen(j,k);c_gen(j,k)],options);
        % [xestimate, functionvalue1]= fmincon(opt_below_diag,[beta_gen(j,k); b_gen(j,k);c_gen(j,k);beta_gen(j,k); 10*b_gen(j,k);c_gen(j,k)],[],[],[],[],[ -1000; 0.1; -1000; -1000; 10; -1000],[ 1000; 80; 2200; 1000; 80; 2200],[],options)
        % For 18 lags
   [xestimate, functionvalue1]= fmincon(opt_below_diag,[beta_gen(j,k); b_gen(j,k);c_gen(j,k);beta_gen(j,k); b_gen(j,k);c_gen(j,k)],[],[],[],[],[ -1000; 20; -1000; -1000; 50; -1000],[ 1000; 30; 5000; 5000; 65; 5000],[],options)

   % Recursive please use fminsearch
    %       [xestimate,functionvalue1]=fminsearch(opt_below_diag,[beta_gen(j,k); b_gen(j,k);c_gen(j,k);beta_gen(j,k); b_gen(j,k);c_gen(j,k)],options);

       
        
        beta_gen12(j,k)=xestimate(1);
        b_gen12(j,k)=xestimate(2);
        c_gen12(j,k)=xestimate(3);
        beta_gen2(j,k)=xestimate(4);
        b_gen2(j,k)=xestimate(5);
        c_gen2(j,k)=xestimate(6);
    end
    
    end
end


%then optimize for three Gaussians

for j=1:setup.size_obs
    for k=1:setup.size_obs
    if k>j %above the diagonal
        opt_restr=@(params)var_quad_restr_3_gaussian( params,setup,store_responses,j,k ) ;
%        [xestimate,functionvalue1]=fminsearch(opt_restr,ones(3,1),options);
       [xestimate,functionvalue1]=fminsearch(opt_restr,[beta_gen(j,k); b_gen(j,k);c_gen(j,k);beta_gen2(j,k); b_gen2(j,k);c_gen2(j,k);beta_gen2(j,k); 2*b_gen2(j,k);2*c_gen2(j,k)],options);
        contemp13(j,k)=0;
        beta_gen13(j,k)=xestimate(1);
        b_gen13(j,k)=xestimate(2);
        c_gen13(j,k)=xestimate(3);
        beta_gen23(j,k)=xestimate(4);
        b_gen23(j,k)=xestimate(5);
        c_gen23(j,k)=xestimate(6);
beta_gen3(j,k)=xestimate(7);
        b_gen3(j,k)=xestimate(8);
        c_gen3(j,k)=xestimate(9);
    elseif j==k
        opt_free=@(params)var_quad_free_2_gaussian( params,setup,store_responses,j,k ) ;
%        [xestimate,functionvalue1]=fminsearch(opt_free,ones(4,1),options);
        [xestimate,functionvalue1]=fminsearch(opt_free,[beta_diag(j);beta_gen(j,k); b_gen(j,k);c_gen(j,k);beta_gen2(j,k); b_gen2(j,k);c_gen2(j,k);beta_gen2(j,k); 2*b_gen2(j,k);2*c_gen2(j,k)],options);
        
        beta_diag13(j)=(xestimate(1));
        
       beta_gen13(j,k)=xestimate(2);
        b_gen13(j,k)=xestimate(3);
        c_gen13(j,k)=xestimate(4);
        beta_gen23(j,k)=xestimate(5);
        b_gen23(j,k)=xestimate(6);
        c_gen23(j,k)=xestimate(7);
        beta_gen3(j,k)=xestimate(8);
        b_gen3(j,k)=xestimate(9);
        c_gen3(j,k)=xestimate(10);
    elseif k<j
    opt_below_diag=@(params)var_quad_below_diag_2_gaussian( params,setup,store_responses,j,k ) ;
%        [xestimate,functionvalue1]=fminsearch(opt_free,ones(4,1),options);
        [xestimate,functionvalue1]=fminsearch(opt_below_diag,[beta_gen(j,k); b_gen(j,k);c_gen(j,k);beta_gen2(j,k); b_gen2(j,k);c_gen2(j,k);beta_gen2(j,k); 2*b_gen2(j,k);2*c_gen2(j,k)],options);
        
       
        
        beta_gen13(j,k)=xestimate(1);
        b_gen13(j,k)=xestimate(2);
        c_gen13(j,k)=xestimate(3);
        beta_gen23(j,k)=xestimate(4);
        b_gen23(j,k)=xestimate(5);
        c_gen23(j,k)=xestimate(6);
       beta_gen3(j,k)=xestimate(7);
        b_gen3(j,k)=xestimate(8);
        c_gen3(j,k)=xestimate(9);
    end
    
    end
end



%I just take the unrestricted initial response as starting value for the
%unrestricted elements of the contemporaneous matrix 
beta_diag=squeeze(setup.store_responses(1:end,setup.index_unrestricted,1));
beta=[];
b=[];
c=[];


%right now only works for up to three Gaussians
beta=[beta;(setup.num_gaussian==1).*beta_gen(:,setup.index_unrestricted)+(setup.num_gaussian==2).*beta_gen12(:,setup.index_unrestricted)+(setup.num_gaussian==3).*beta_gen13(:,setup.index_unrestricted)];

b=[b;(setup.num_gaussian==1).*b_gen(:,setup.index_unrestricted)+(setup.num_gaussian==2).*b_gen12(:,setup.index_unrestricted)+(setup.num_gaussian==3).*b_gen13(:,setup.index_unrestricted)];
%c=[c;(setup.num_gaussian==1).*c_gen(:,setup.index_unrestricted)+(setup.num_gaussian>1).*c_gen12(:,setup.index_unrestricted)];
%c=[c;c_gen2((setup.num_gaussian>1),setup.index_unrestricted)];       
c=[c;(setup.num_gaussian==1).*c_gen(:,setup.index_unrestricted)+(setup.num_gaussian==2).*c_gen12(:,setup.index_unrestricted)+(setup.num_gaussian==3).*c_gen13(:,setup.index_unrestricted)];
% if max(setup.num_gaussian)>1
beta=[beta;(setup.num_gaussian==2).*beta_gen2(:,setup.index_unrestricted)+(setup.num_gaussian==3).*beta_gen23(:,setup.index_unrestricted)];
beta=[beta;(setup.num_gaussian==3).*beta_gen3(:,setup.index_unrestricted)];
b=[b;(setup.num_gaussian==2).*b_gen2(:,setup.index_unrestricted)+(setup.num_gaussian==3).*b_gen23(:,setup.index_unrestricted)];
b=[b;(setup.num_gaussian==3).*b_gen3(:,setup.index_unrestricted)]; 
c=[c;(setup.num_gaussian==2).*c_gen2(:,setup.index_unrestricted)+(setup.num_gaussian==3).*c_gen23(:,setup.index_unrestricted)];
c=[c;(setup.num_gaussian==3).*c_gen3(:,setup.index_unrestricted)]; 
  
% end 

 
beta2=beta(beta~=0);
b2=b(b~=0);
c2=c(c~=0);
beta_temp=reshape(beta,setup.size_obs,3);
beta_temp=sum(beta_temp,2);
 
temp_length=length(c);


%setup.initial_parameter=[zeros(setup.size_obs,1);beta_diag;beta2(:);b2(:);c2(:);beta_diag;beta2(:);b2(:);c2(:)];
setup.initial_parameter=[mean(data,2);beta_diag;beta2(:);b2(:);c2(:);beta_diag;beta2(:);b2(:);c2(:)];




%plot fitted and VAR IRFs
epsilon_vec=zeros(setup.size_obs,setup.lags+1);
[ Sigma, intercept] = unwrap_NL_IRF( setup.initial_parameter,epsilon_vec,setup,setup.indicator );
%for now I only plot the response to the unrestricted shock
figure;
for jj=1:setup.size_obs
   subplot(setup.size_obs,1,jj)
   plot(1:size(Sigma,3),squeeze(Sigma(jj,setup.index_unrestricted,:)),1:size(Sigma,3),squeeze(setup.store_responses(jj,setup.index_unrestricted,:)))
end



