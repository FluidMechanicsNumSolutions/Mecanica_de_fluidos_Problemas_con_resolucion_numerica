function [ M2 ] = MNumNormShock( M1, gamma )
%Dado un número de Mach (supersonico) calcula el número de Mach 
%detrás de la onda de choque plana
% M1: Numero de Mach aguas arriba de la onda de choque
% M2: Numero de Mach aguas abajo de la onda de choque
% gamma: Índice de politropía
M2=((M1^2+2/(gamma-1))/(2*gamma*M1^2/(gamma-1)-1))^(1/2);
end

