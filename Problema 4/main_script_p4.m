%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCRIPT PRINCIPAL DEL PROBLEMA 4 (v1.0.0)

%%%% Los resultados que se obtienen en este programa
%%%% son los asociados al problema 4 del libro:
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
close all
clc

% Definición de la geometría, propiedades físicas y condiciones de contorno
R=0.5; % Radio R (m)
avect=[0.07 0.085 0.1]; % Radio a (m) ==>Se incluyen los diferentes radios 
b=0.05; % Radio b (m)
S=pi*b^2; % Superficie S (m^2)
den1=8000; % Densidad del material del conducto (kg/m^3)
rho=1000; % Densidad del fluido (kg/m^3)
Qo=0.0035; % Caudal inicial (m^3/s)
k=10; % Constante k (N·m·s)

i=1; 
for a=avect
    fprintf('Grosor de pared: %.2f cm (a=%.3f m)\n',(a-b)*100,a)
    % Cálculo de la velocidad angular en régimen permanente 
    wrp = 4*rho*R*S*((Qo*exp(5)/(2*S))^2)/(k+8*rho*S*R^2*Qo*exp(5)/(2*S)); 
    fprintf('\t - Velocidad angular en reg. permanente: %.2f rad/s\n',wrp)
    % Cálculo de la masa del conducto y del momento de inercia
    mcond=2*den1*pi^2*R*(a^2-b^2); 
    Iz=(1/4)*(8*R^2+3*a^2+3*b^2)*mcond; 
    fprintf('\t - Inercia en eje z: %.2f kg*m^2\n',Iz)
    % Resolución de la ecuación diferencial
    param=[R a b S den1 mcond rho Qo k Iz]; 
          % Parámetros de entrada en ode45
    % -->Condicion inicial para la primera ecuación diferencial
    wo=0;
    CI=wo; 
    % -->Tiempo final (Régimen permanente) 
    tfinal=5;


    % -->Solución de la primera ecuación diferencial
    [t1,wrt1]=ode45( @(t,x) f1(t,x,param),linspace(0,tfinal,1000),CI);
    % --->Condicion inicial para la segunda ecuación diferencial
    CI=wrt1(end);
    to=tfinal;
    % -->Solución de la segunda ecuación diferencial
    [t2,wrt2]=ode45( @(t2,x) f2(t2,x,param,to),...
        linspace(tfinal,tfinal+5,1000),CI);
    
    % Se unen los vectores para obtener la solución total
    t=[t1; t2];
    w=[wrt1; wrt2];
    

    % Cálculo del tiempo para llegar a régimen permanente 
    [ti,index]=unique(t);
    trp=interp1(w(index),ti,wrp*0.998);
    fprintf(['\t - Tiempo reg. perm. con error ' ...
        'del 0.2%%: %.2f segundos\n'],trp)
    % Guarda datos
    save(['Datos_caso_' num2str(i)])
    % Representación gráfica
    hold on;
    
    if a==0.070
        figure(1)
        p{1}=plot(t1,wrt1,t2,wrt2,'LineWidth',1.5);
        ylabel('w (rad/s)');
        xlabel('t (s)');
        legend('ec. dif. 1','ec. dif. 2')
        grid minor
        ax=gca;
        ax.FontSize=12;
        p{1}(1).LineStyle='--';
        p{1}(2).LineStyle='-.';

    end
     figure(2)
     hold on; p{2}(i)=plot(t,w,'LineWidth',1.5);
     ylabel('w (rad/s)');
     xlabel('t (s)');
     ax=gca;
     ax.FontSize=12;
   
i=i+1;
end
figure (2); legend(['a=' num2str(avect(1)) ' m'],...
    ['a=' num2str(avect(2)) ' m'],['a=' num2str(avect(3)) ' m'])
p{2}(1).Marker='^'; p{2}(1).MarkerIndices=(25:250:length(t));  
p{2}(2).Marker='*'; p{2}(2).MarkerIndices=(185:250:length(t));
p{2}(3).Marker='o'; p{2}(3).MarkerIndices=(75:250:length(t));
grid minor

