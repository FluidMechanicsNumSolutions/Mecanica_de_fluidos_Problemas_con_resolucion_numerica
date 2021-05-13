function AP = distr_P_down(Gamma,rho,X,Xo,Yo)
%Funci�n para obtener la distribuci�n de presi�n sobre la placa

% Uinf: Velocidad no perturbada
% (Xo, Yo): Posici�n del v�rtice
% Gamma: Circulaci�n asociada al v�rtice
% X: Coordenada horizontal de la posici�n a calcular
% rho: Densidad del aire
% AP: Presi�n relativa

    A=X.^2.*(-Xo^2 + Yo^2) - (Xo^2 + Yo^2).^2;
    B=((X - Xo).^2 + Yo^2).*((X + Xo).^2 + Yo^2);
    AP=-rho*Gamma^2/(2*pi^2).*((4*X*Xo*Yo./B).^2 + A./(B.*(Xo^2 + Yo^2)));
    
end

