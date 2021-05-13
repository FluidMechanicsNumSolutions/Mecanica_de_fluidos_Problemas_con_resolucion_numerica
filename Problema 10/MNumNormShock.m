function [ M2 ] = MNumNormShock( M1, gamma )
%Dado un numero de Mach (supersonico) calcula el numero de Mach detras de la onda de choque 
% M1: Numero de Mach aguas arriba de la onda de choque
% M2: Numero de Mach aguas abajo de la onda de choque
% gamma: Índice de politropía
M2=((M1^2+2/(gamma-1))/(2*gamma*M1^2/(gamma-1)-1))^(1/2);
end

