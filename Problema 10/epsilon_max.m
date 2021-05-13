function epsilon_maxima = epsilon_max(M1,gamma)
% Función que calcula el ángulo épsilon máximo 
%   M1: Número de Mach entrante
%   gamma: Constante politrópica
%   epsilon_maxima: Ángulo épsilon máximo

epsilon_maxima=asin(sqrt(((gamma+1)/4*M1^2-1+((gamma+1)*(1+...
    (gamma-1)/2*M1^2+(gamma+1)/16*M1^4))^(1/2))/(gamma*M1^2)));
end


