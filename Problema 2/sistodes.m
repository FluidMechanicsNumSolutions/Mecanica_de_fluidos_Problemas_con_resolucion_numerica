function [ dxdt ] = sistodes(t,x,Dem)
% Derivadas de primer orden 
% t: tiempo
% x: Variables que se integrar�n: 
%    => x(1): velocidad del �mbolo
%    => x(2): posici�n del �mbolo 
%    => x(3): Presi�n de la c�mara de aceite
%    => x(4): Presi�n de la c�mara de nitr�geno
% dxdt: Derivadas de primer orden de cada una de las variables: 
%    => dxdt(1): aceleraci�n del �mbolo
%    => dxdt(2): velocidad del �mbolo 
%    => dxdt(3): variaci�n de la presi�n de la c�mara de aceite respecto del
%                tiempo
%    => dxdt(4): variaci�n de la presi�n de la c�mara de nitr�geno respecto del
%                tiempo

global Sem Dp Sp ho Vp B gamma m g VN2in PN2in  Pacin 

u=x(1);
y=x(2);
Pac=x(3);
PN2=x(4);
Vp=-(ho/2)*sin(t);
dxdt(1)= (-PN2*Sem+Pac*Sem-m*g)/m;
dxdt(2)=u;
dxdt(3)=(-Sp*Vp-Sem*u)*B/(Sem*y+(1/2)*Sp*ho*(1+cos(t)));
dxdt(4)=Sem*u*gamma*(PN2^((1/gamma)+1))/(VN2in*PN2in^(1/gamma));

dxdt=dxdt(:);
end
