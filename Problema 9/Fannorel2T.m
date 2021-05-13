function [ relT1T2 ] = Fannorel2T( M,gamma )
% En flujo de Fanno, relaciona la temperatura estatica entre dos puntos
% M: Número de Mach en el punto 1
% gamma: Índice de politropía
% relT1T2=T1/T2 donde el flujo va de 1 a 2 (Atencion: no al reves)
% En el punto 2 hay condicion critica M2=1;
relT1T2=((gamma+1)*(1/2))/(1+((gamma-1)*(1/2))*M^2);
end

