function f= f_eta(gamma,M,choose_M_eta,eta_inf)
% Función necesaria para obtener el ángulo eta o el Número de Mach
% en una onda de expansión de Prandtl-Meyer
%   gamma: Índice de politropía
%   M: Número de Mach entrante 
%   f: Variable de salida
%   choose_M_eta: 
%           => Si choose_M_eta es "eta"
%              La variable de salida f será eta (ángulo de Prandtl-Meyer)
%           => Si choose_M_eta es "M"
%              La variable de salida f será el Número de Mach "M"
%   eta_inf: Ángulo eta inferior partiendo de un Numero de Mach unitario   
 

eta=sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)...
        *((M.^2-1))./(gamma+1)))-atan(sqrt((M.^2-1))); %Ángulo de Prandtl-Meyer
if choose_M_eta=="eta"
    f=eta;
else
    if choose_M_eta=="M"
       AM=1e-5;
       Msolve= (1:AM:100);
       f=sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)...
        *((Msolve.^2-1))./(gamma+1)))-atan(sqrt((Msolve.^2-1)))-eta_inf;
       for i=1:length(Msolve)-1
             if f(i)*f(i+1)<0 && Msolve(i)>M
                sol=(Msolve(i)+Msolve(i+1))/2; 
             end
       end
       f=sol;
    end
end

