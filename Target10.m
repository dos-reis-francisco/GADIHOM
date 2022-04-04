%% data values for Target 7
nchromosomes=128; % number of chromosomes (must be 4 factor)
nkeep=64;   % number of keeped chromosomes
rhov=0.1; % homogenized volumic density target
seed=6; % number of beams per side 
nkmax=1000; % number iterations max
% Material properties target
Ex=3000; Ey=3000; Gxy=1000; etayxy=-1.0; etaxxy=-1.0; nuyx=-1.0;
target=Compliance(Ex, Ey, Gxy, etayxy, etaxxy, nuyx); % compliance tensor
wtarget=[1 1 1 10 10 10]; %weight vector
mutrate=0.05;  % mutation rate
ntvalue=1000;   % number of different value for beam width t
lambda=0.0;   % weight of  
convergence=0.001;
nConvergence=50;
GADIHOM(nchromosomes,nkeep,rhov,seed,nkmax,target,wtarget,mutrate,ntvalue,...
    lambda,convergence,nConvergence);