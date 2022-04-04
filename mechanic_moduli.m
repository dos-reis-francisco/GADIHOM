function [K,Ex,Ey,nuyx,nuxy,muxy,etaxxy,etayxy,etaxyx,etaxyy] = mechanic_moduli(MS)
K=1/(MS(1,1)+MS(1,2)+MS(1,6)+MS(1,6));
Ex=1/MS(1,1);
Ey=1/MS(1,2);
nuyx=-MS(1,6)*Ey;
nuxy=-MS(1,6)*Ex;
muxy=1/(2*MS(1,3));
etaxxy=MS(1,5)*2*muxy;
etayxy=MS(1,4)*2*muxy;
etaxyx=MS(1,5)*Ex;
etaxyy=MS(1,4)*Ey;
end

