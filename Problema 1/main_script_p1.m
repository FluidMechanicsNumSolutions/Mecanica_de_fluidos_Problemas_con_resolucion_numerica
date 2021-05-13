
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCRIPT PRINCIPAL DEL PROBLEMA 1 (v1.0.0)

%%%% Los resultados que se obtienen en este programa
%%%% son los asociados al problema 1 del libro:
%%%%
%%%% MECÁNICA DE FLUIDOS. PROBLEMAS CON RESOLUCIÓN NUMÉRICA.
%%%% Autores: Álvaro Mañas González y Josep M. Bergada Granyó 
%%%%
%%%% Para mayor información acerca del enunciado del problema 
%%%% y la resolución del mismo, es necesario adquirir el libro completo

%%%% Este script es el principal, el que debe ser ejecutado
%%%% para obtener los resultados definidos en el problema.
%%%% Los datos de entrada pueden ser modificados. 
%%%% Una vez ejecutado este script, se presentan los resultados 
%%%% en Command Window y en las figuras corresponientes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
close all

global Rp Di Sp alpha Ao Vo l1 l2 l4 l6 l8 l10 l3 l5 l7 l9...
       l11 hs1 hs3 hs2 h1 h10 h2 h11 r1 r2 r3 r4...
       Pd w Cd B mu rho Ptank
 
% Geometría
% => La nomenclatura utilizada es la misma que en problema
% => Las unidades utilizadas están en Sistema Internacional

Rp=0.060; 
D=0.014;
Di=D;
Sp=pi*((Di/2)^2);
Spi=pi*((0.008/2)^2);
alpha=20*pi/180;
Do=0.006;
Ao=pi*((Do/2)^2);
Vo=Spi*0.04+Sp*0.0005;
l1=0.002;
l2=0.0007;
l4=l2;
l6=l2;
l8=l2;
l10=l2;
l3=0.004;
l5=l3;
l7=l3;
l9=l3;
l11=0.019;
hs1=10*10^(-6);
hs3=hs1;
hs2=(0.4*10^(-3))+hs1;
h1=2*10^(-6); 
h10=0.0004+h1; h2=h10; h3=h1; h4=h10; h5=h1; h6=h1;
h7=h10; h8=h1; h9=h10; h11=h1;
r1=0.005;
r2=0.0075;
r3=0.008;
r4=0.010;

% Definción de las condiciones de contorno y propiedades físicas
% => La nomenclatura utilizada es la misma que en problema
% => Las unidades utilizadas están en Sistema Internacional
Pd=300*10^5;
Cd=0.6;
B=1.8*10^9;
tin=0;
mu=857*32*10^(-6);
pin=300*10^5;
rho=857;
Ptank=10^5;
w=2000*2*pi/60;
tfinal=pi/w;

%Resolución numérica 

%%%%%% PISTÓN ASCENDENTE 
%   - Con fugas y sin fugas
%   - Caso inicial con hs=10 micras y h1=2 micras

[t,pres_conf,QL1,QL2,t1s,pres_sinf] = Res_num_subida(tfinal,pin);
save('sin_con_QLh10')


%   - Con fugas
%   - Con otras distancias entre slipper y plato 
%      [ hs=20 micras y hs=25 micras ]
%   - Se mantiene h1=2 micras

for hs1=[20,25]*10^(-6)
    hs3=hs1;
    hs2=(0.4*10^(-3))+hs1;
    [t,pres_conf,QL1,QL2] = Res_num_subida(tfinal,pin);
    save(['con_QLh' num2str(hs1*10^(6))])
end

%   - Con fugas
%   - Con otras distancias entre cilindro y piston
%     [ h=3 micras , h=4 micras y h=5 micras ]
%   - Se mantiene hs=10 micras

hs1=10*10^(-6);
hs3=hs1;
hs2=(0.4*10^(-3))+hs1;

for h1=[3,4,5]*10^(-6)
    h10=0.0004+h1; h2=h10; h3=h1; h4=h10; h5=h1; 
    h6=h1; h7=h10; h8=h1; h9=h10; h11=h1;
    [t,pres_conf,QL1,QL2] = Res_num_subida(tfinal,pin);
    save(['con_H' num2str(h1*10^(6))])
end

%%%%%% PISTÓN DESCENDENTE

hs1=10*10^(-6);
hs3=hs1;
hs2=(0.4*10^(-3))+hs1;
h1=2*10^(-6);
h10=0.0004+h1; h2=h10; h3=h1; h4=h10; h5=h1; h6=h1; 
h7=h10; h8=h1; h9=h10; h11=h1;

%----> Si las condiciones son: 
Pd=1*10^5;
pin=1*10^5;
[t,p_down, QL1] = Res_num_bajada(tfinal,pin);
save('down_1bar')

%----> Si se aumenta la presión en la entrada para evitar la cavitación: 
Pd=7.5*10^5;
pin=7.5*10^5;
[t,p_down, QL2] = Res_num_bajada(tfinal,pin);
save('down_7_5bar')

%Representación gráfica 
plots; % Ejecución del script 'plots.m'
