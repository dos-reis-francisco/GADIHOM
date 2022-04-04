%% Compliance.m
% Compute the compliance tensor from mechanics moduli
% Dos Reis F.
% 17.02.2021
function target= Compliance(Ex, Ey, Gxy, etayxy, etaxxy, nuyx )
S11=1/Ex;
S22=1/Ey;
S66=1/(2*Gxy);
S26=etayxy/(2*Gxy);
S16=etaxxy/(2*Gxy);
S12=-nuyx/Ey;
target=[S11, S22, S66, S26, S16, S12] ;
end

