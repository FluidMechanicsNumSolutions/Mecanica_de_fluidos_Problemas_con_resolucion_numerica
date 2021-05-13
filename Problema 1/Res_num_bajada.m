function [t,p_down,QL1,QL2] = Res_num_bajada(tfinal,pin)
% Resolución numérica del sistema de ecuaciones en el descenso
% pin: Presión inicial
% tfinal: Tiempo final
% t: vector de tiempo
% p_down: vector de presión durante la bajada
% QL1: caudal volumétrico circulante de fuga entre pistón y cilindro
% QL2: caudal volumétrico circulante de fuga entre slipper y plato
% inclinado

global Rp Di Sp alpha Ao Vo l1 l2 l4 l6 l8 l10 l3 l5 l7 l9...
       l11 hs1 hs3 hs2 h1 h10 h2 h11 r1 r2 r3 r4...
       Pd w Cd B mu rho Ptank
%-->CÁLCULO DE PRESIONES
%Resolución de la ecuacion diferencial: 

[t,p_down]=ode23('f3_down',linspace(0,tfinal,1000),pin);
p_down=real(p_down);
%-->CÁLCULO DE CAUDALES DE FUGAS
Q=zeros(1,length(t));
QL1=zeros(1,length(t));
QL2=zeros(1,length(t));

for i=1:length(t)
    p=p_down(i);
    %p=u3(i);
    tc=t(i);
    QL1(i)=pi*Di*((h1*Rp*tan(alpha)*sin(w*tc)*w)*(1/2))-...
        (1/12)*pi*Di*(Ptank-p-(6*Rp*tan(alpha)*(sin(w*tc)))*...
        w*mu*(h10-h1)*(l2+l4+l6+l8+l10)/h10^3)/(mu*(l1+l2+l3+...
        l4+l5+l6+l7+l8+l9+l10+l11-(0.195*10^(-1))*(1/2)-...
        Rp*tan(alpha)*cos(w*tc))/(h11^3)+mu*(1/h2^3-1/h1^3)*...
        (l2+l4+l6+l8+l10));
    QL2(i)=pi*(p-Ptank)/((6*mu)*(log(r2/r1)/hs1^3+...
        log(r3/r2)/hs2^3+log(r4/r3)/hs3^3));

end
end

