function [Mflow] =massflow(M,Po,To,S,R,gamma)
%Flujo m�sico en un punto dado el M,Po,To y gamma 
%   M: Mach number
%   Po: Presi�n de estancamiento [Pa]
%   To: Temperatura de estancamiento [K]
%   gamma: �ndice de politrop�a
%   Mflow: Flujo m�sico [kg/s]

Mflow=((gamma/(R*To))^0.5)*Po*S*M*(1+((M^2)*(gamma-1)/2))^...
    (-(gamma+1)/(2*(gamma-1)));

end

