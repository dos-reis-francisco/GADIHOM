%% fcost.m
% cost function for the genetic algorithm 
% sub module for inverse homogenization code
% Dos Reis F.
% 02.2021
function cost = fcost(A,B,w,lambda)
normeA=sqrt(A*A');
normeB=sqrt(B*B');
D1=A/normeA;
D2=B/normeB;
D=D1-D2;
D3=D.*w;
cost=sqrt(D3*D3')+lambda*normeB/normeA;
end

