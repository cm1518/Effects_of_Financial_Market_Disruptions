% Main prg to run the BM asymmetric estimation

clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part to change depending on whether RB or CM is running the code
fileroot='C:/Users/Alexander Ziegenbein/Dropbox/asymmetry_EBP';
addpath([fileroot,'/dataandlibrary'])  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialization/Parametrization of routine
setup_NL;
%set indicator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting values
%To find starting values, match IRFs with those of VAR:
   VAR_resp_match_lp;
close all;

setup.length_param_vector=length(setup.initial_parameter);
%scaling for adaptive MCMC (see handbook of MCMC, page 104) ADAPTED
setup.scaling_adaptive=.03^2/setup.length_param_vector;

setup.index_normal=[5; 6; 33; 34];
setup.index_gamma=[];
%% Recursive identification via tight priors: no effect of negative EBP shocks (positive or negative) on GDP and PI 
setup.normal_prior_means=[zeros(4,1)];
setup.normal_prior_std=[ones(4,1)*0.001];

%estimation
[ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH( setup );

%add_matrices are the estimated shocks
save results