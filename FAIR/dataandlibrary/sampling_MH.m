function [ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH( setup )
%main function for the MH algorithm
data=load(setup.data);
data=1*data.data;

%line below has been added!!
data=data(:,end-setup.sample_size+1:end);
%



warning off;
%options = optimset('Display','iter','TolFun',1e-12,'TolX',1e-8,'MaxIter',1e5, 'MaxFunEvals',1e6);
options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-4,'MaxIter',100, 'MaxFunEvals',50000);

if setup.proposal==1
   temp=setup.scaling_draws/setup.check_scaling;
   %pr_function=str2func(setup.prior_function);
   
   if setup.likelihood==1
      % ll_function=SSKF_wrap( params,setup,data );
       post=@(params,setup,data) (prior(params,setup)+SSKF_wrap( params,setup,data ));
       postopt=@(params,setup,data) -(prior(params,setup)+SSKF_wrap( params,setup,data ));
       postopt2=@(params) -(prior(params,setup)+SSKF_wrap( params,setup,data ));
       %since we use a minimizer
   elseif setup.likelihood==2
       %ll_function=KF_wrap;
        post=@(params,setup,data) (prior(params,setup)+KF_wrap( params,setup,data ));
       postopt=@(params,setup,data) -(prior(params,setup)+KF_wrap( params,setup,data )); %since we use a minimizer
       postopt2=@(params) -(prior(params,setup)+KF_wrap( params,setup,data ));
   elseif setup.likelihood==3
    post=@(params,setup,data) (prior(params,setup)+likelihood_wrap( params,setup,data ));
       postopt=@(params,setup,data) -(prior(params,setup)+likelihood_wrap( params,setup,data )); %since we use a minimizer
       postopt2=@(params) -(prior(params,setup)+likelihood_wrap( params,setup,data ));
   end
%    post=@(params) (pr_function(params,setup)+ll_function(params,setup,data));
%    postopt=@(params) -(pr_function(params,setup)+ll_function(params,setup,data)); %since we use a minimizer
if setup.skip_opt==0
   %initial minimization
   [ x0 ] = transform( setup.initial_parameter,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub  );
 
   
   [xestimate1,functionvalue1]=fminsearch(postopt2,x0,options);
    options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-4,'MaxIter',300, 'MaxFunEvals',50000);

  [xestimate1,functionvalue1,exitflag,output,grad,hessian]=fminunc(postopt2,xestimate1,options);
   
   %Hesstemp = hessian2(func2str(postopt),xestimate,setup,data);
   %Hesstemp=reshape(Hesstemp,setup.length_param_vector,setup.length_param_vector);
  % [functionvalue3, xestimate3,gh,H,itct,fcount,retcodeh] =csminwel(postopt,xestimate2,hessian,[],1e-10,50,setup,data);
  % if functionvalue1==min([functionvalue1,functionvalue2,functionvalue3])
       xestimate=xestimate1;
       functionvalue=functionvalue1;
  % elseif functionvalue2==min([functionvalue1,functionvalue2,functionvalue3])
  % xestimate=xestimate2;
  % functionvalue=functionvalue2;
  % else
  %     xestimate=xestimate3;
  %     functionvalue=functionvalue3;
  % end
 save temp_max xestimate1
    elseif setup.skip_opt==1
load temp_max
xestimate=xestimate1;
functionvalue=postopt2(xestimate)
end
  
  
%   
%   Hess = hessian2(func2str(post),xestimate,setup,data);
%    %fdhess(func2str(post),xestimate,setup,data)
%    Hess=reshape(-Hess,setup.length_param_vector,setup.length_param_vector);
%    Hesstemp=eye(size(Hess,1))/Hess;
%    
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%    %making sure the Hessian is positive definite
%   try 
%    [V,D]=eig(Hesstemp);
% 
%    d=diag(D);
%    d(d<=0)=eps;
% 
%    Hess= V*diag(d)*V';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   catch ME
      %disp('Could not compute Hessian, will use scaled identity matrix instead.')
      
%I will set the initial values for the covariances by hand here
 for bb=1:setup.number_blocks      
Hess{bb}=.5*eye(length(setup.index_block{bb}));
end
  %end
   scaling=ones(setup.number_blocks,1).*setup.initial_scaling;
   old_posterior=-functionvalue;
   posterior_draw=xestimate;
   acc_rate=zeros(setup.number_blocks,1);
   for jj=1:temp %loop to adjust scaling matrix
       disp('sample for scaling matrix')
       jj
    for bb=1:setup.number_blocks
       [ acc_rate(bb) ,posterior_draw,old_posterior] = acceptance_rate( posterior_draw,Hess{bb},scaling(bb),setup.check_scaling,old_posterior,setup,data,bb  );  
     acc_rate(bb) 
     if acc_rate(bb)<.2
       scaling(bb)=.75*scaling(bb);
       elseif acc_rate(bb)>.5
           scaling(bb)=1.5*scaling(bb);
     end
    end 
    end
    disp('scaling:')
    scaling
    disp('acceptance rate:')
    acc_rate
    
    
%     %burn_in
%        [ acc_rate ,posterior_draw,old_posterior] = acceptance_rate( posterior_draw,Hess,scaling,setup.burn_in,old_posterior,setup,data  ); 
%     
       
    %actual draws
    
    [ acc_rate ,draws,log_posteriors, statedraws, add_matrices] = draws_standard_RW( posterior_draw,Hess,scaling,old_posterior,setup,data  );
    
 else %adaptive MH
     
     disp('this code is only implemented for the standard MH algorithm for now')
%   
%      temp=setup.scaling_draws/setup.check_scaling;
%    %pr_function=str2func(setup.prior_function);
%    
%    if setup.likelihood==1
%       % ll_function=SSKF_wrap( params,setup,data );
%        post=@(params,setup,data) (prior(params,setup)+SSKF_wrap( params,setup,data ));
%        postopt=@(params,setup,data) -(prior(params,setup)+SSKF_wrap( params,setup,data ));
%        postopt2=@(params) -(prior(params,setup)+SSKF_wrap( params,setup,data ));
%        %since we use a minimizer
%    elseif setup.likelihood==2
%        %ll_function=KF_wrap;
%         post=@(params,setup,data) (prior(params,setup)+KF_wrap( params,setup,data ));
%        postopt=@(params,setup,data) -(prior(params,setup)+KF_wrap( params,setup,data )); %since we use a minimizer
%        postopt2=@(params) -(prior(params,setup)+KF_wrap( params,setup,data ));
%     elseif setup.likelihood==3
%     post=@(params,setup,data) (prior(params,setup)+likelihood_wrap( params,setup,data ));
%        postopt=@(params,setup,data) -(prior(params,setup)+likelihood_wrap( params,setup,data )); %since we use a minimizer
%        postopt2=@(params) -(prior(params,setup)+likelihood_wrap( params,setup,data ));
%    end
% %    post=@(params) (pr_function(params,setup)+ll_function(params,setup,data));
% %    postopt=@(params) -(pr_function(params,setup)+ll_function(params,setup,data)); %since we use a minimizer
%    %initial minimization
%   if setup.skip_opt==0
%    
%    [ x0 ] = transform( setup.initial_parameter,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub  );
%    [xestimate1,functionvalue1]=fminsearch(postopt2,x0,options);
%    options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-4,'MaxIter',100, 'MaxFunEvals',4000); 
%     %RB: try the fminunc (sometimes crash the code (with bad initial guesses(?)) so use a try)
%     try
%         disp('running fminunc')
%         [xestimate2,functionvalue2,exitflag,output,grad,hessian]=fminunc(postopt2,xestimate1,options);
%         %Hesstemp = hessian2(func2str(postopt),xestimate,setup,data);
%         %Hesstemp=reshape(Hesstemp,setup.length_param_vector,setup.length_param_vector);
%         %[functionvalue3, xestimate3,gh,H,itct,fcount,retcodeh] =csminwel(postopt,xestimate2,hessian,[],1e-10,50,setup,data);
%         
%         if functionvalue1==min([functionvalue1,functionvalue2])
%             xestimate=xestimate1;
%             functionvalue=functionvalue1;
%         elseif functionvalue2==min([functionvalue1,functionvalue2])
%             xestimate=xestimate2;
%             functionvalue=functionvalue2;
%             %else
%             %    xestimate=xestimate3;
%             %    functionvalue=functionvalue3;
%         end
%     catch me
%         xestimate=xestimate1;
%         functionvalue=functionvalue1;
%         disp('  ')
%         disp('fminunc could not be used')
%         disp('  ')
%           [functionvalue3, xestimate3,gh,H,itct,fcount,retcodeh] =csminwel(postopt,xestimate,hessian,[],1e-10,50,setup,data);
%         if functionvalue3<functionvalue1
%             xestimate=xestimate3;
%         elseif functionvalue3>functionvalue1
%             xestimate=xestimate1;
%         end
%     end
%     
%   elseif setup.skip_opt==1
% xestimate=setup.initial_parameter;
%    [ xestimate ] = transform( xestimate,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub  );
%    functionvalue=postopt2(xestimate)
%   end
% %    disp('Calculating Hessian')
% %    Hess = hessian2(func2str(post),xestimate,setup,data);
% %    Hess=reshape(-Hess,setup.length_param_vector,setup.length_param_vector);
% %    Hesstemp=eye(size(Hess,1))/Hess;
% %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %    %making sure the Hessian is positive definite
% % 
% %   try 
% %    [V,D]=eig(Hesstemp);
% % 
% %    d=diag(D);
% %    d(d<=0)=eps;
% % 
% %    Hess= V*diag(d)*V';
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   catch ME
%       disp('Could not compute Hessian, will use scaled identity matrix instead.')
%      Hess=0.00001*diag(abs(xestimate)); 
%      
%      save interim xestimate functionvalue
%   %end
%   % disp('Done calculating Hessian')
%    scaling=setup.initial_scaling;
%    old_posterior=-functionvalue;
%    posterior_draw=xestimate;
%    for jj=1:temp %loop to adjust scaling matrix
%        disp('sample for scaling matrix')
%        jj
%      [ acc_rate ,posterior_draw,old_posterior] = acceptance_rate( posterior_draw,Hess,scaling,setup.check_scaling,old_posterior,setup,data  );  
%      acc_rate 
%      if acc_rate<.2
%        scaling=.75*scaling;
%        elseif acc_rate>.5
%            scaling=1.5*scaling;
%        end
%    end
%     disp('scaling:')
%     scaling
%     disp('acceptance rate:')
%     acc_rate
%     
%     
%  
%         
%        
%     %burn_in
%     
%      [ acc_rate ,draws,log_posteriors, statedraws] = draws_standard_RW_for_Adaptive( posterior_draw,Hess,scaling,old_posterior,setup,data  );
%     
%      meandraws=mean(draws,2);
%      cov_initial=(setup.scaling_adaptive^(-1))*scaling*Hess-setup.eps_adaptive*eye(setup.length_param_vector);
%       %first_draw=transform( draws(:,end),setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub  );
%      %cov_initial=0.00000000000000000000001^2*abs(diag(meandraws));
%     %actual draws
%     
%     [ acc_rate ,draws,log_posteriors,statedraws, add_matrices] = draws_adaptive_RW( draws(:,end),cov_initial,log_posteriors(end),setup,data ,meandraws );
%     
end



end

