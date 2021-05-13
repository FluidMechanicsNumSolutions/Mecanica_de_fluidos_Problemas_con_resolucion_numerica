function [ rPoP ] = relPoP( M, gamma )
%Rel. entre presi�n de estancamiento con presi�n est�tica EN UN PUNTO  
%M: N�mero de mach en el punto 
%gamma: �ndice de politrop�a
%rPoP: Relaci�n entre la presi�n de estancamiento y 
       %la presi�n est�tica en el punto considerado

rPoP=(1+(1/2)*M^2*gamma-(1/2)*M^2)^(gamma/(gamma-1));

end

