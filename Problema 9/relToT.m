function [ rToT ] = relToT( M,gamma )
%Rel. entre temperatura de estancamiento con temperatura est�tica EN UN PUNTO  
%M: N�mero de mach en el punto 
%gamma: �ndice de politrop�a
%rToT: Relaci�n entre la temperatura de estancamiento y
       %la temperatura est�tica en el punto considerado

rToT=(M^2*gamma-M^2+2)/2;
end

