%% data values for Target 2
nchromosomes=128; % number of chromosomes (must be 4 factor)
nkeep=64;   % number of keeped chromosomes
rhov=0.1; % homogenized volumic density target
seed=6; % number of beams per side 
nkmax=1000; % number iterations max
% Material properties target
Ex=3000; Ey=3000; Gxy=1500; etayxy=0.0; etaxxy=0.0; nuyx=-0.3;
target=Compliance(Ex, Ey, Gxy, etayxy, etaxxy, nuyx); % compliance tensor
wtarget=[10 10 10 10 10 50]; %weight vector
mutrate=0.05;  % mutation rate
ntvalue=1000;   % number of different value for beam width t
lambda=4.0;   % weight of  
convergence=0.003;
nConvergence=50;
GADIHOM(nchromosomes,nkeep,rhov,seed,nkmax,target,wtarget,mutrate,ntvalue,...
    lambda,convergence,nConvergence);