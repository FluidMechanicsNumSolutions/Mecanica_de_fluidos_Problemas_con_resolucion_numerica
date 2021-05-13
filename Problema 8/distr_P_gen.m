function AP = distr_P_gen(Uinf,Gamma,rho,X,Y,Xo,Yo,V_g)
%Funci�n para obtener la distribuci�n de presi�n en un punto gen�rico

% Uinf: Velocidad no perturbada
% (Xo, Yo): Posici�n del v�rtice
% Gamma: Circulaci�n asociada al v�rtice
% X,Y: Coordenadas de la posici�n a calcular
% rho: Densidad del aire
% V_g: Velocidad en m�dulo en el punto gen�rico 
% AP: Presi�n relativa

dphidx=Gamma./(2*pi).*(-(Y+Yo)./((X-Xo).^2+(Y+Yo).^2)+...
    (Y-Yo)./((X-Xo).^2+(Y-Yo).^2)+(Y-Yo)./((X+Xo).^2+...
    (Y-Yo).^2)-(Y+Yo)./((X+Xo).^2+(Y+Yo).^2));

dphidy=Gamma./(2*pi).*(-(X-Xo)./((X-Xo).^2+(Y+Yo).^2)-...
    (X-Xo)./((X-Xo).^2+(Y-Yo).^2)+(X+Xo)./((X+Xo).^2+...
    (Y-Yo).^2)+(X+Xo)./((X+Xo).^2+(Y+Yo).^2));

U=Uinf+(Gamma/(4*pi))*(Xo^2)/(Yo*(Xo^2 + Yo^2));
V=-(Gamma/(4*pi))*Yo^2/(Xo*(Xo^2 + Yo^2));

dphidt=dphidx*U+dphidy*V;

AP=-rho*(dphidt+V_g.^2/2);

end

