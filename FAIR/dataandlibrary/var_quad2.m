function [ VAR_match ] = var_quad2( params,setup,store_responses,j,k )
%objective function for VAR responses matching
[ MA_matrices ] = wrap_BM_var_opt_11( params, setup );
MA_h=MA_matrices(j,k,1:11);
sr=store_responses(j,k,:);
% size(MA_h)
% size(store_responses)

VAR_match=(MA_h(:)-sr(:))'*(MA_h(:)-sr(:));
end

