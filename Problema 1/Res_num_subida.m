function [t,p_conf,QL1,QL2,ts,p_sinf] = Res_num_subida(tfinal,Pin)
% Resolución numérica del sistema de ecuaciones en el ascenso
% pin: Presión inicial
% tfinal: Tiempo final
% t: vector de tiempo (Caso con fugas)
% p_conf: vector de presión durante la subida (Caso con fugas)
% QL1: caudal volumétrico circulante de fuga entre pistón y cilindro
% QL2: caudal volumétrico circulante de fuga entre slipper y plato
% inclinado
% ts: vector de tiempo (Caso sin fugas)
% p_sinf: vector de presión durante la subida (Caso sin fugas)

global Rp Di Sp alpha Ao Vo l1 l2 l4 l6 l8 l10 l3 l5 l7 l9...
       l11 hs1 hs3 hs2 h1 h10 h2 h11 r1 r2 r3 r4...
       Pd w Cd B mu rho Ptank

%----Sin fugas
%-->CÁLCULO DE PRESIONES
%Resolución de la ecuacion diferencial: 
[t2,u2]=ode23('f2',linspace(0,0.0007833,1000),Pin);
[t21,u21]=ode23s('f2',linspace(0.0007833,0.01428,1000),u2(length(u2)));
[t22,u22]=ode23('f2',linspace(0.01428,tfinal,1000),u21(length(u21)));
ts=[t2' t21' t22']';
p_sinf=real([u2' u21' u22']');


%----Con fugas
%-->CÁLCULO DE PRESIONES
%Resolución de la ecuación diferencial:  
[t4,u3]=ode23('f3_up',linspace(0,0.0008725,1000),300*10^5);
[t5,u4]=ode23s('f3_up',linspace(0.0008725,0.01428,1000),u3(length(u3)));
[t6,u5]=ode23('f3_up',linspace(0.01428,tfinal,1000),u4(length(u4)));
t=[t4' t5' t6']';
p_conf=real([u3' u4' u5']');

%-->CÁLCULO DE CAUDALES DE FUGAS
Q=zeros(1,length(t));
QL1=zeros(1,length(t));
QL2=zeros(1,length(t));

for i=1:length(t)
    p=p_conf(i);  
    tc=t(i);
    Q(i)=Cd*Ao*sqrt(2*(p-Pd)/rho);
    QL1(i) = pi*Di*((-h1*Rp*tan(alpha)*sin(w*tc)*w)*(1/2))-...
        (1/12)*pi*Di*(Ptank-p-(6*Rp*tan(alpha)*(-sin(w*tc)))...
        *w*mu*(h10-h1)*(l2+l4+l6+l8+l10)/h10^3)/(mu*...
        (l1+l2+l3+l4+l5+l6+l7+l8+l9+l10+l11-(0.195*10^(-1))*...
        (1/2)-Rp*tan(alpha)*cos(w*tc))/(h11^3)+mu*...
        (1/h2^3-1/h1^3)*(l2+l4+l6+l8+l10));
    QL2(i) = pi*(p-Ptank)/((6*mu)*((log(r2/r1)/(hs1^3))+...
        (log(r3/r2)/(hs2^3))+(log(r4/r3)/(hs3^3))));
end


end

