function [ relMach ] = FannoIndvar( M,gamma )
% Calcula 4CfL/Dh conociendo el numero de Mach .
% Recuerda que 4Cf=frictionfactor=f
% M: Número de Mach en el punto 1
% gamma: Índice de politropía
relMach=(-M^2+1)/(gamma*M^2)+(gamma+1)*log((gamma+1)*...
    M^2/(2*(1+((gamma-1)*(1/2))*M^2)))/(2*gamma);
end

