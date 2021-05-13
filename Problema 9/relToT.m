function [ rToT ] = relToT( M,gamma )
%Rel. entre temperatura de estancamiento con temperatura estática EN UN PUNTO  
%M: Número de mach en el punto 
%gamma: Índice de politropía
%rToT: Relación entre la temperatura de estancamiento y
       %la temperatura estática en el punto considerado

rToT=(M^2*gamma-M^2+2)/2;
end

