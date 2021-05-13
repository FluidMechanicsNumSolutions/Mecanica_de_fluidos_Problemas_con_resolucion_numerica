function [Mflow] =massflow(M,Po,To,S,R,gamma)
%Flujo másico en un punto dado el M,Po,To y gamma 
%   M: Mach number
%   Po: Presión de estancamiento [Pa]
%   To: Temperatura de estancamiento [K]
%   gamma: Índice de politropía
%   Mflow: Flujo másico [kg/s]

Mflow=((gamma/(R*To))^0.5)*Po*S*M*(1+((M^2)*(gamma-1)/2))^...
    (-(gamma+1)/(2*(gamma-1)));

end

