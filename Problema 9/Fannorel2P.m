function [ relP1P2 ] = Fannorel2P( M,gamma )
% En flujo de Fanno, relaciona la presion estatica entre dos puntos
% M: Número de Mach en el punto 1
% gamma: Índice de politropía
% relP1P2=P1/P2 donde el flujo va de 1 a 2 (Atencion no al reves)
% En el punto 2 hay condicion critica M2=1;
relP1P2=(((gamma+1)*(1/2))/(1+((gamma-1)*(1/2))*M^2))^(1/2)/M;
end

