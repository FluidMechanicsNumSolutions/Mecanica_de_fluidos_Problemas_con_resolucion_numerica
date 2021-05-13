function [ relPo1Po2 ] = Fannorel2Po( M,gamma )
% En flujo de Fanno, relaciona la presion de estancamiento entre dos puntos
% M: N�mero de Mach en el punto 1
% gamma: �ndice de politrop�a
% relPo1Po2=Po1/Po2 donde el flujo va de 1 a 2 (Atenci�n: no al rev�s)
% En el punto 2 hay condici�n cr�tica M2=1;
relPo1Po2=((1+((gamma-1)*(1/2))*M^2)/((gamma+1)*...
    (1/2)))^((gamma+1)/(2*(gamma-1)))/M;
end