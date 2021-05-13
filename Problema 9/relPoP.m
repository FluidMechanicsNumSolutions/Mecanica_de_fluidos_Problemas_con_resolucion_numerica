function [ rPoP ] = relPoP( M, gamma )
%Rel. entre presión de estancamiento con presión estática EN UN PUNTO  
%M: Número de mach en el punto 
%gamma: Índice de politropía
%rPoP: Relación entre la presión de estancamiento y 
       %la presión estática en el punto considerado

rPoP=(1+(1/2)*M^2*gamma-(1/2)*M^2)^(gamma/(gamma-1));

end

