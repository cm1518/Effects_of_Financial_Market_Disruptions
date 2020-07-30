%creates cell with objects needed for estimation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%number of draws 
setup.number_of_draws=100000;% in the paper, we use 2,000,000
%length of burn in (cafter calibrating scaling factor for the case of
%standard RW proposal)
%setup.burn_in=1000;%10000
%number of draws for choosing the scaling matrix in case of standard RW
%proposalfigur
setup.scaling_draws=4000;%10000
%scaling is recomputed after check_scaling draws
setup.check_scaling=100;%200
%initial scaling for standard RW proposal
setup.initial_scaling=[.00001 0.5 50 1]';
setup.proposal=1;
setup.likelihood=3;
%name of function that evaluates log prior
setup.prior_function='prior'; %right now prior function has to be called prior.m!
%initial value for the state and covariance of the state
setup.skip_opt=0; %skip optimization and go directly to MCMC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setup.lags=120; %lag length of IRF, 10 years
setup.size_obs=4; %number of observables
setup.freq=1;% frequency of data: 1 for monthly, 2 for quarterly
setup.shocks=0; %0-> initial shocks are zero, 1-> initial shocks from VAR (reduces sample size)
setup.polynomials=0; %degree of polynomial detrending
setup.symmetry=1000; %making sure non of the options are turned on
setup.VARsym_order=12; %order of symmetric VAR used for starting values
%that's the number6f lags in the sym VAR (RB)

%data file
run (['DGPSR',num2str(setup.size_obs),'_BM'])
run (['symVAR_estimation']);
save ((['VAR',num2str(setup.size_obs),'.mat']), 'A', 'sigma','errors')
setup.data=['data_file',num2str(setup.size_obs),'.mat'];

%setting up initial shocks - either form VAR or set to zero
if setup.shocks==0
  setup.initial_eps=zeros(setup.size_obs,setup.lags); %intial shocks set to 0 
  setup.sample_size=length(data);
else
    setup.initial_eps=errors(:,1:setup.lags);
    setup.sample_size=length(data)-setup.VARsym_order-setup.lags;
end

%display every disp iter draw
setup.disp_iter=1000;
% keep every keep_draw th draw
setup.keep_draw=10;

%size of state vector
setup.state_size=1;
%epsilon for adaptive RW
%setup.eps_adaptive=.005^2;
setup.eps_adaptive=.5^2;

%initial parameter value for optimization
load ((['VAR',num2str(setup.size_obs),'.mat']))

setup.VARsymA=A; % A matrix for option 3 above
setup.VARsymcov=sigma; %covariance matrix of residuals from estimated VAR (for option 3 above)
% If recursive VAR;
%setup.VARsymchol=chol(sigma,'lower'); 

% If proxy SVAR
setup.VARsymchol = A_full; % 

%priors
setup.index_normal=[];
setup.index_gamma=[];

%setting up restrictions 
setup.index_restricted=[1 2 4];
setup.index_unrestricted=[ 3 ];

%% Some restrictions to keep the parameters from "running away" and speed up convergence
% I restrict b to be smaller than 60months or 5 years, i.e.
% the peak effects of a shock have to occur within the first 5 years.
% I use the lower bounds to assure that the Gaussian function don't lie on
% top of each other, but rather explore the IR at differenz horizons.
% The bounds avoid that some parameters drift into the "wrong" direction (low
% likelihood) and thus speed up convergence.

setup.index_logit_general=[17:24 45:52];
setup.length_logit_general=length(setup.index_logit_general);
setup.logit_general_lb=[-15*ones(1,3) 10 -15*ones(1,3) -10 -15*ones(1,8) ]';
setup.logit_general_ub=[60*ones(1,16)]';

setup.length_log=2;
setup.index_log=[7 35]; % all 2GB

setup.length_logit=0;
setup.index_logit=[];

setup.threshold_vec=1.65*ones(3,1);

%should additional matrices be stored
setup.add_matrices=1;
%dimension of those variables
setup.dim_add_matrices=[setup.size_obs setup.sample_size+setup.lags];
setup.number_blocks=4;
setup.index_block{1}=[1:4 7:8 35:36 37:44 9:16]'; % constants, contemporaneous impact, a
setup.index_block{2}=[17:24 45:52]'; % b
setup.index_block{3}=[25:32 53:60]'; % c
setup.index_block{4}=[5:6 33:34]'; % restricted parameters (see Principal.m)

setup.num_gaussian=[2 2 2 2]'; % 2 Gaussians for each IR

%% Redundant if we do not allow for state-dependence
load indicator
setup.indicator=indicator;