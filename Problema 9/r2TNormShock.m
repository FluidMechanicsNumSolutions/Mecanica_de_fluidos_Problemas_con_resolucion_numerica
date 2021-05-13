function [ r2T ] = r2TNormShock( M1, gamma )
% Relaci�n de temperaturas est�ticas aguas abajo y arriba de la onda
% Conociendo el n�mero de Mach del fluido situado aguas arriba de la onda,
% mediante esta funci�n se obtiene la relaci�n de temperaturas: T2/T1
% M1: N�mero de Mach aguas arriba de la onda de choque
% gamma: �ndice de politrop�a
% r2T: Relaci�n de temperaturas est�ticas (T2/T1)
r2T=(1+((gamma-1)*(1/2))*M1^2)*(2*gamma*M1^2/(gamma-1)-1)*(2*(gamma-1))/...
    ((gamma+1)^2*M1^2);
end

