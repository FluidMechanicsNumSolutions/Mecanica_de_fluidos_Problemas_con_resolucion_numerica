%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCRIPT PRINCIPAL DEL PROBLEMA 6 (v1.0.0)

%%%% Los resultados que se obtienen en este programa
%%%% son los asociados al problema 6 del libro:
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

%%%%% Apartados 1 y 2
% Geometría y propiedades físicas
% => La nomenclatura utilizada es la misma que en problema
% => Las unidades utilizadas están en Sistema Internacional

Rr=0.01;
R=0.01;
Ro=Rr+4*10^(-6);
H=4*10^(-6);
L=3*10^(-3);
gamma=40*pi/180;
delta=45*pi/180;
Rp=Rr+(R+H)*sin(gamma)-(R+H)*sin(gamma+delta);
mu=1.983*10^(-5);
params=[Rr R H gamma];

% Condiciones de contorno
% => La nomenclatura utilizada es la misma que en problema
% => Las unidades utilizadas están en Sistema Internacional

Pin=100*10^5;
Pout=1*10^5;

% Resolución de la integral Idelta
Idelta=integral(@(alpha) dI(alpha,params),0,delta);
% Constantes halladas 
A=(log(Ro/Rr)*Ro^4-log(Ro/Rr)*Rr^4-Ro^4+2*Ro^2*Rr^2-Rr^4)/log(Ro/Rr);
P1=(-4*L*H^3*Pout/((3*A*(2*R+H))*Idelta)+Pin)/...
    (1-4*L*H^3/((3*A*(2*R+H))*Idelta));
Qcon=pi*H^3*(Pout-P1)/(6*mu*(2*R+H)*Idelta);
Qcil=-A*pi*(P1-Pin)/(8*mu*L);
% Definición de la distribución de presión en las dos zonas
%--->Zona cilíndrica
x=(0:0.001:L);
Pcil=(P1-Pin).*x/L+Pin;
%--->Zona cónica
ang=(0:0.0001:delta);
Pcon=zeros(1,length(ang));
Pcon(1)=P1;
for i=2:length(ang)
    Ialpha=integral(@(alpha) dI(alpha,params),0,ang(i)...
        ,'RelTol',1e-10,'AbsTol',1e-10);
    Pcon(i)=(Pout-P1)*Ialpha/Idelta+P1;
end  

% Fuerza horizontal debida a la presión
Fx=0; %Definición de la variable
for i=1:length(ang)-1
    % Este for se utiliza para realizar una integración numérica
    % Discretización:
    angprom=(ang(i)+ang(i+1))/2; 
    Pconprom=(Pcon(i)+Pcon(i+1))/2; 
    difang=ang(i+1)-ang(i);
    % Fuerza realizada sobre cada superficie discretizada
    Fxcalc=Pconprom*cos(gamma+angprom)*2*pi*(R+H)*...
        (Rr+(R+H)*(sin(gamma)-sin(gamma+angprom)))*difang;
    % Suma de fuerzas entre todas las superficies:
    Fx=Fx+Fxcalc;% Fuerza de presión
end
fprintf('Caudal zona cilíndrica: %.3e m^3/s \n',Qcil)
fprintf('Caudal zona cónica: %.3e m^3/s \n',Qcon)
fprintf('Fuerza de presión: %.3f N \n',Fx);

% Representación gráfica
figure(1)
figure(gcf)
plot(x,Pcil,'LineWidth',1.2)
xlabel('Coordenada horizontal (m)')
ylabel('Presión (Pa). Zona cilíndrica.')
ax=gca;
ax.FontSize=12;
figure(2)
figure(gcf)
ang=ang*180/pi;
plot(ang,Pcon,'LineWidth',1.2)
xlabel('Ángulo \alpha (^o)')
ylabel('Presión (Pa). Zona cónica.')
ax=gca;
ax.FontSize=12;

%%%%% Apartado 3

% Geometría

%R: Hace referencia al radio del cilindro, pero
%el valor sigue siendo el mismo
h=10*1e-6; 
e=1e-6; 


% Fuerza horizontal y vertical sobre el eje
delta=linspace(0,2*pi,1000);
w=5000; 
[Fx_T,Fy_T] = Fx_T_Fy_T(delta,e,mu,w,R,h);

figure(3)
plot(delta*180/pi,Fx_T,'LineWidth',1.2); 
xlabel('Ángulo \delta (º)'); 
ylabel('Fuerza horizontal F_{X´T} [N/m]'); 
ylim([-3000 3000]); xlim([0 360]);
ax=gca;
ax.FontSize=12;ax.XTick=(0:60:360);

figure(4)
plot(delta*180/pi,Fy_T,'LineWidth',1.2); 
xlabel('Ángulo \delta (º)'); 
ylabel('Fuerza vertical F_{Y´T} [N/m]'); 
ylim([-3000 3000]); xlim([0 360]);
ax=gca;
ax.FontSize=12; ax.XTick=(0:60:360);  

% Par restaurador total
ang_casos=[45, 225; ...
           90, 270; ...
          135, 315]; %Casos de estudio
      
e_max=e; 

% Discretización del eje cilíndrico

n=500; % Número de porciones
i=1:n; % Indices
e_vect=e_max*(1-(i-1)/n-1/2/n); %Excentricidad en el centro
                                %de cada porción
d=1/2-(i-1)/2/n-1/n/4; %Distancia respecto el centro
e_cent=[e_vect rot90(e_vect,2)]; %Vector total de excentricidades
d=[d rot90(d,2)]; %Vector total de distancias
Lx=linspace(-1/2+1/n/2,1/2-1/n/2,2*n); %Coordenada longitudinal
ang_cent=zeros(1,2*n);% Ángulo delta en el centro de cada porción

% Inicialización de variables
w_vect=linspace(0,5000,11); %Velocidades angulares
Mz_tot=zeros(3,length(w_vect)); 
plots=cell(1,3);

% Cálculo de los tres casos estudiados
linstyle=["-","-","--"];
for caso=1:3

    ang_extrem=[ang_casos(caso,1),ang_casos(caso,2)];
    ang_cent(1:n)=ang_extrem(1)*pi/180;
    ang_cent(n+1:end)=ang_extrem(2)*pi/180;

    Fyi=zeros(1,n);
    Mzi=zeros(1,n);

  for j=1:length(w_vect) 
    w=w_vect(j);
    for i=1:length(d)
        [Fx_T,Fy_T] = Fx_T_Fy_T(ang_cent(i),e_cent(i),mu,w,R,h);
        Fyi(i)=Fy_T*(d(1)-d(2));
        if ang_cent(i)>pi/2 && ang_cent(i)<3*pi/2
           Mzi(i)=-Fyi(i)*d(i);
        else
           Mzi(i)=+Fyi(i)*d(i);
        end
    end
    Mz_tot(i,j)=sum(Mzi);
    
    if w_vect(j)==5000
        
        figure(5)
        plot(Lx,Fyi,'LineWidth',1.2); hold on; 

        figure(6)
        plot(Lx,Mzi,linstyle(caso),'LineWidth',1.2); hold on;
        
    end
  end
  
  figure(7)
  plots{caso}=plot(w_vect,Mz_tot(i,:),'-o'); hold on;
  
end

figure(5)
legend('Caso 45^o-225^o','Caso 90^o-270^o',...
            'Caso 135^o-315^o')
xlabel('Coordenada longitudinal al eje (m)');
ylabel('Fuerza F_{Yi} (N/m)')
title('Fuerza F_{Yi} sobre cada porción (\omega = 5000 rad/s)')
ax=gca;
ax.FontSize=12;

figure(6)
legend('Caso 45^o-225^o','Caso 90^o-270^o',...
            'Caso 135^o-315^o')
xlabel('Coordenada longitudinal al eje (m)');
ylabel('Momento generado por F_{Yi} (N·m)');
title('Momento M_{Zi} sobre cada porción (\omega = 5000 rad/s)')
ax=gca;
ax.FontSize=12;

figure(7)
plots{1}.MarkerIndices=(1:2:length(w_vect));
plots{3}.MarkerIndices=(2:2:length(w_vect));
legend('Caso 45^o-225^o','Caso 90^o-270^o',...
            'Caso 135^o-315^o')
xlabel('Velocidad angular \omega (rad/s)');
ylabel('Par restaurador M_{ZoT} (N·m)/m');
ax=gca;
ax.FontSize=12;
