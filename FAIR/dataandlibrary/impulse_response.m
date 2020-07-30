function [ IRFS ] = impulse_response(shock, sign, draw,setup )
%returns the impoulse response to a certain shock, the sign of which also
%has to be specified:shock is the index of the shock and sign is 0 if the
%shock is negative and 1 if it is positive

ind=sign+1;

IRFS=zeros(setup.size_obs,setup.lags+1);

[ constant current_matrices lagged_matrices ] = wrap_BM_irf( draw, setup );

IRFS(:,1)=current_matrices(:,shock,ind);

for jj=1:setup.lags
  IRFS(:,jj+1)= lagged_matrices(:,shock,ind,jj); 
    
end
end

