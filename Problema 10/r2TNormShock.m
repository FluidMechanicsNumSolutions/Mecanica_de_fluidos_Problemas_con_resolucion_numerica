function [ r2T ] = r2TNormShock( M1, gamma )
% Relación de temperaturas estáticas aguas abajo y arriba de la onda
% Conociendo el número de Mach del fluido situado aguas arriba de la onda,
% mediante esta función se obtiene la relación de temperaturas: T2/T1
% M1: Número de Mach aguas arriba de la onda de choque
% gamma: Índice de politropía
% r2T: Relación de temperaturas estáticas (T2/T1)
r2T=(1+((gamma-1)*(1/2))*M1^2)*(2*gamma*M1^2/(gamma-1)-1)*(2*(gamma-1))/...
    ((gamma+1)^2*M1^2);
end

