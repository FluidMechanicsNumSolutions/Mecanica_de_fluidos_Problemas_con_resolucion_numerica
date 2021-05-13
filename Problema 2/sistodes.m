function [ dxdt ] = sistodes(t,x,Dem)
% Derivadas de primer orden 
% t: tiempo
% x: Variables que se integrarán: 
%    => x(1): velocidad del émbolo
%    => x(2): posición del émbolo 
%    => x(3): Presión de la cámara de aceite
%    => x(4): Presión de la cámara de nitrógeno
% dxdt: Derivadas de primer orden de cada una de las variables: 
%    => dxdt(1): aceleración del émbolo
%    => dxdt(2): velocidad del émbolo 
%    => dxdt(3): variación de la presión de la cámara de aceite respecto del
%                tiempo
%    => dxdt(4): variación de la presión de la cámara de nitrógeno respecto del
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
