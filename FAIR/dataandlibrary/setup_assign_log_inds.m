function [ setup ] = setup_assign_log_inds( setup )
%assigns the log indices in the setup structure (tedious to do by hand)

setup.length_log=2*setup.size_obs;

ind=zeros(setup.size_obs+length(setup.index_asymmetric_vars),1);

%matrix for negative shocks

ind(1)=setup.size_obs+1;


for jj=1:(setup.size_obs-1)
ind(jj+1)=ind(jj)+setup.size_obs-(jj-1);
end


%matrix for positive shocks

for jj=1:length(setup.index_asymmetric_vars)
    if jj==1
    ind(jj+setup.size_obs)=ind(jj+setup.size_obs-1)+1;
    elseif jj>1
ind(jj+setup.size_obs)=ind(jj+setup.size_obs-1)+1+(setup.size_obs-setup.index_asymmetric_vars(jj-1));
    end
end




setup.index_log=ind;


end

