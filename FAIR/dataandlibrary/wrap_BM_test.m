function [ constant current_matrices lagged_matrices ] = wrap_BM_test( params, setup )
%function that maps a parameter vector into the objects needed to evaluate
%the likelihood function

%transforming parameters back to the original parameter space
[ params ] = inv_transform( params,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );

constant=params(1:setup.size_obs);

current_matrices=zeros(setup.size_obs,setup.size_obs,2);

current_matrices(:,:,1)=ltvec(params(setup.size_obs+1:setup.size_obs+(setup.size_obs*(setup.size_obs+1))/2));
current_matrices(:,:,2)=ltvec(params(setup.size_obs+(setup.size_obs*(setup.size_obs+1))/2+1:setup.size_obs+setup.size_obs*(setup.size_obs+1)/2+(setup.size_obs*(setup.size_obs+1))/2));

lagged_matrices=zeros(setup.size_obs,setup.size_obs,2,setup.lags);


ind_for_loop=setup.size_obs+setup.size_obs*(setup.size_obs+1)/2+(setup.size_obs*(setup.size_obs+1))/2;
for jj=1:setup.lags_no_symmetry
   lagged_matrices(:,:,1,jj)=reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs); 
   ind_for_loop=ind_for_loop+setup.size_obs^2;
   lagged_matrices(:,:,2,jj)=reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs); 
   ind_for_loop=ind_for_loop+setup.size_obs^2;
end


if setup.symmetry==1 %no conditions beyond symmetry imposed
   for kk=1:(setup.lags-setup.lags_no_symmetry) 
    lagged_matrices(:,:,1,kk+setup.lags_no_symmetry)=reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs);
    lagged_matrices(:,:,2,kk+setup.lags_no_symmetry)=reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs);
    ind_for_loop=ind_for_loop+setup.size_obs^2;
   end

elseif setup.symmetry==2 %VAR(1) type symmetry imposed
    A=reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs);
    ind_for_loop=ind_for_loop+setup.size_obs^2;
    Sigma=ltvec(params(ind_for_loop+1:end));
    for kk=1:(setup.lags-setup.lags_no_symmetry) 
      lagged_matrices(:,:,1,kk+setup.lags_no_symmetry)=(A)^(kk+setup.lags_no_symmetry)*Sigma;
      lagged_matrices(:,:,2,kk+setup.lags_no_symmetry)=lagged_matrices(:,:,1,kk+setup.lags_no_symmetry);
    end
    
elseif setup.symmetry==3 %estimates from VAR(1) imposed
for kk=1:(setup.lags-setup.lags_no_symmetry) 
 lagged_matrices(:,:,1,kk+setup.lags_no_symmetry)=(setup.VARsymA)^(kk+setup.lags_no_symmetry)*setup.VARsymchol;
 lagged_matrices(:,:,2,kk+setup.lags_no_symmetry)=lagged_matrices(:,:,1,kk+setup.lags_no_symmetry);
end

elseif setup.symmetry==4 %VAR imposed on coefficient matrices
    phi=zeros(setup.size_obs^2,setup.size_obs^2,setup.order_VAR_A); %get coefficient matrices for the VAR law of motion of the original coefficient matrices
    
    for vv=1:setup.order_VAR_A
       phi(:,:,vv)=reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^4),setup.size_obs^2,setup.size_obs^2);
       ind_for_loop=ind_for_loop+setup.size_obs^4;
    end
    
    
   for kk=1:(setup.lags-setup.lags_no_symmetry) 
 
    if setup.lags_no_symmetry+1==setup.order_VAR_A %first lag that enters the law of motion for the coefficents is the contemporaneous response
    
    temp=setup.implied_VAR_matrices(:,kk+setup.order_VAR_A);
    for aa=1:setup.order_VAR_A
        if (kk==1) && (aa == setup.order_VAR_A)
        temp2=reshape(current_matrices(:,:,1),setup.size_obs^2,1);
        temp=temp+phi(:,:,aa)*(temp2-setup.implied_VAR_matrices(:,kk+setup.order_VAR_A-aa));
        
        else
        temp2=reshape(lagged_matrices(:,:,1,kk+setup.lags_no_symmetry-aa),setup.size_obs^2,1);
        temp=temp+phi(:,:,aa)*(temp2-setup.implied_VAR_matrices(:,kk+setup.order_VAR_A-aa));
        end    
    end
    
    
    
  lagged_matrices(:,:,1,kk+setup.lags_no_symmetry)=reshape(temp,setup.size_obs,setup.size_obs);
 
 
 temp=setup.implied_VAR_matrices(:,kk+setup.order_VAR_A);
    for aa=1:setup.order_VAR_A
        if (kk==1) && (aa == setup.order_VAR_A)
        temp2=reshape(current_matrices(:,:,2),setup.size_obs^2,1);
        temp=temp+phi(:,:,aa)*(temp2-setup.implied_VAR_matrices(:,kk+setup.order_VAR_A-aa));
        
        else
        temp2=reshape(lagged_matrices(:,:,2,kk+setup.lags_no_symmetry-aa),setup.size_obs^2,1);
        temp=temp+phi(:,:,aa)*(temp2-setup.implied_VAR_matrices(:,kk+setup.order_VAR_A-aa));
        end    
    end
    
    
    
  lagged_matrices(:,:,2,kk+setup.lags_no_symmetry)=reshape(temp,setup.size_obs,setup.size_obs);
 
 
 
 
 
 
    
 
    elseif setup.lags_no_symmetry+1>setup.order_VAR_A %contemporaneous response matrix never enters law of motion for coefficients
        temp=setup.implied_VAR_matrices(:,kk+setup.order_VAR_A);
        disp('test')
    for aa=1:setup.order_VAR_A
        temp2=reshape(lagged_matrices(:,:,1,kk+setup.lags_no_symmetry-aa),setup.size_obs^2,1);
        temp=temp+phi(:,:,aa)*(temp2-setup.implied_VAR_matrices(:,kk+setup.order_VAR_A-aa));
    end
    
    
    lagged_matrices(:,:,1,kk+setup.lags_no_symmetry)=reshape(temp,setup.size_obs,setup.size_obs);
 
    temp=setup.implied_VAR_matrices(:,kk+setup.order_VAR_A);
    for aa=1:setup.order_VAR_A
        temp2=reshape(lagged_matrices(:,:,2,kk+setup.lags_no_symmetry-aa),setup.size_obs^2,1);
        temp=temp+phi(:,:,aa)*(temp2-setup.implied_VAR_matrices(:,kk+setup.order_VAR_A-aa));
    end
 
 
 
 
    lagged_matrices(:,:,2,kk+setup.lags_no_symmetry)=reshape(temp,setup.size_obs,setup.size_obs);

    
    end
    end


end  
    
end

