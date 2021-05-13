function [ relPo1Po2 ] = Fannorel2Po( M,gamma )
% En flujo de Fanno, relaciona la presion de estancamiento entre dos puntos
% M: Número de Mach en el punto 1
% gamma: Índice de politropía
% relPo1Po2=Po1/Po2 donde el flujo va de 1 a 2 (Atención: no al revés)
% En el punto 2 hay condición crítica M2=1;
relPo1Po2=((1+((gamma-1)*(1/2))*M^2)/((gamma+1)*...
    (1/2)))^((gamma+1)/(2*(gamma-1)))/M;
end