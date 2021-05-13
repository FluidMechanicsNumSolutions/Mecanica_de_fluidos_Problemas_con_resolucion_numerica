function [ M2 ] = MNumNormShock( M1, gamma )
%Dado un n�mero de Mach (supersonico) calcula el n�mero de Mach 
%detr�s de la onda de choque plana
% M1: Numero de Mach aguas arriba de la onda de choque
% M2: Numero de Mach aguas abajo de la onda de choque
% gamma: �ndice de politrop�a
M2=((M1^2+2/(gamma-1))/(2*gamma*M1^2/(gamma-1)-1))^(1/2);
end

