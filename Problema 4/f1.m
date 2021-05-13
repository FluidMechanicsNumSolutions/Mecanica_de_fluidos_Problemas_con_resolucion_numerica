function [ dxdt ] = f1( t,x,param )
% Derivadas de primer orden
% t: tiempo
% x: Variables que se integrarán: 
%    => x(1): velocidad angular
% param: parámetros necesarios definidos en el script principal
% dxdt: Derivadas de primer orden de cada una de las variables
%    => dxdt(1): aceleración angular

R=param(1);a=param(2);b=param(3);S=param(4);den1=param(5);m1=param(6);
rho=param(7);Qo=param(8);k=param(9); Iz= param(10);
w=x(1);
dxdt(1)=(-(-4*rho*R^2*Qo*exp(t)-k)*w-4*rho*R*S*(Qo*exp(t)/(2*S))^2-...
    rho*pi*R^2*Qo*exp(t))/(-4*pi*R^3*S*rho-Iz);
dxdt=dxdt(:);
end

