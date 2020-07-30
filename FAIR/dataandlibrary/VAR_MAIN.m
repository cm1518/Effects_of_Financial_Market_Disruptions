clear all
clc


load ..\Data\data_file3.mat %data has to be stored in an array called data.

%calling Koop's BVAR code (using independent Normal_Wishart.)
p = 24;               % Number of lags on dependent variables

BVAR_GIBBS;

%reshape parameter draws

ALPHA_draws_reshaped =reshape(ALPHA_draws,nsave,K*M); 
for jj=1:nsave;
SIGMA_draws_reshaped(jj,:)= unique(SIGMA_draws(jj,:,:));
end

draws=[ALPHA_draws_reshaped SIGMA_draws_reshaped];

md_var=m_harmonic(draws,log_posteriors);
marg_likelihood=exp( m_harmonic(draws,log_posteriors));