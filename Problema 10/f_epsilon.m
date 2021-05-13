function sol = f_epsilon(delta,gamma,M1)
% Función necesaria para obtener las raices de la ecuación que permite
% determinar el ángulo épsilon. 
%   Epsilon: Ángulo de la onda de choque oblicua respecto el eje
%               horizontal
%   delta: Inclinación del sólido
%   gamma: Índice de politropía
%   M1: Número de Mach entrante 
%   sol: Raiz hallada de la ecuación (Solución Epsilon)



Aepsilon=1e-5; % Incremento del ángulo

Epsilon_min=asin(1/M1); % Ángulo mínimo sin presencia de discontinuidad 
Epsilon_max=epsilon_max(M1,gamma); % Ángulo máximo
Epsilon=(Epsilon_min+Aepsilon:Aepsilon:Epsilon_max);


f = ((gamma+1)./2.*M1.^2./(M1.^2.*sin(Epsilon).^2-1)-1).*...
    tan(Epsilon)-1./tan(delta);
sol=zeros(1,2);
k=0;
for i=1:length(Epsilon)-1
   
    if f(i)*f(i+1)<0 && Epsilon(i)>delta 
       k=k+1;
       sol(k)=(Epsilon(i)+Epsilon(i+1))/2; 
    end
    
end
if k<2
    sol=sol(1);
end
end

