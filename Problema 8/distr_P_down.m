function AP = distr_P_down(Gamma,rho,X,Xo,Yo)
%Función para obtener la distribución de presión sobre la placa

% Uinf: Velocidad no perturbada
% (Xo, Yo): Posición del vórtice
% Gamma: Circulación asociada al vórtice
% X: Coordenada horizontal de la posición a calcular
% rho: Densidad del aire
% AP: Presión relativa

    A=X.^2.*(-Xo^2 + Yo^2) - (Xo^2 + Yo^2).^2;
    B=((X - Xo).^2 + Yo^2).*((X + Xo).^2 + Yo^2);
    AP=-rho*Gamma^2/(2*pi^2).*((4*X*Xo*Yo./B).^2 + A./(B.*(Xo^2 + Yo^2)));
    
end

