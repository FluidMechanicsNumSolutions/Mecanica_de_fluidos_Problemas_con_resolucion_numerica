%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCRIPT PRINCIPAL DEL PROBLEMA 2 (v1.0.0)

%%%% Los resultados que se obtienen en este programa
%%%% son los asociados al problema 2 del libro:
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
clc
clear all

global Sem Dp Sp ho Vp B gamma m g VN2in PN2in  Pacin 

%Tres casos
Dem1=130*10^(-3);
Dem2=120*10^(-3);
Dem3=110*10^(-3);
Dems=[Dem1 Dem2 Dem3];

%Propiedades físicas
% => La nomenclatura utilizada es la misma que en problema
% => Las unidades utilizadas están en Sistema Internacional

B=1.5*10^9;
gamma=1.4;
m=5;
g=9.81;
for i=1:3
    indice=i;
    %Geometría
    Dem=Dems(i);
    Sem=pi*(Dem^2)/4;
    Dp=0.055;
    Sp=pi*(Dp^2)/4;
    ho=0.2;
    %Condiciones de contorno y Condiciones iniciales
    % => La nomenclatura utilizada es la misma que en problema
    % => Las unidades utilizadas están en Sistema Internacional

    VN2in=Sem*0.1;
    PN2in=10*10^5;
    Pacin=(PN2in*Sem+m*g)/Sem;
    w=1;
    tfinal=2*pi/w;
    yemin=0.65;
    uin=0;
    CI=[uin,yemin,Pacin,PN2in];
    % Resolución del sistema de ecuaciones diferenciales
    [t,y]=ode23s(@(t,x) sistodes(t,x,Dem),[0 2*tfinal], CI);
    %Guarda los resultados 
    if i==1
        save('Datos_1'); 
    elseif i==2
        save('Datos_2');
    else
        save('Datos_3');
    end
end
% Representación gráfica
plots % Ejecución del script 'plots.m' 
