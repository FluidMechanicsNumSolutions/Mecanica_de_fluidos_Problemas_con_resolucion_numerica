function [ r2P ] = r2PNormShock( M1, gamma )
% Relacion de presiones estaticas aguas abajo y arriba de la onda
% Conociendo el numero de Mach del fluido situado aguas arriba de la onda,
% mediante esta funcion se obtiene la relacion de presiones estaticas: P2/P1 
% M1: Numero de Mach aguas arriba de la onda de choque
% gamma: Índice de politropía
% r2P: Relacion de presiones estáticas (P2/P1)

r2P=2*gamma*M1^2/(gamma+1)-(gamma-1)/(gamma+1);
end
