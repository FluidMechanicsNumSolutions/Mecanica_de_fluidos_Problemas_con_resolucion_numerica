function epsilon_maxima = epsilon_max(M1,gamma)
% Funci�n que calcula el �ngulo �psilon m�ximo 
%   M1: N�mero de Mach entrante
%   gamma: Constante politr�pica
%   epsilon_maxima: �ngulo �psilon m�ximo

epsilon_maxima=asin(sqrt(((gamma+1)/4*M1^2-1+((gamma+1)*(1+...
    (gamma-1)/2*M1^2+(gamma+1)/16*M1^4))^(1/2))/(gamma*M1^2)));
end


