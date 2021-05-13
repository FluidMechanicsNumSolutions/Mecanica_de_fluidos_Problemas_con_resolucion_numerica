close all
clear all

% Previamente, es necesario tener los datos de las integraciones
% guardados si SOLO se desea ver las gráficas, teniendo los
% cálculos hechos. Si se ejecuta el archivo principal, 
% este script se ejecuta de forma automática. 

fsize=11;
%Representación gráfica (CASO:PISTÓN ASCENDENTE)
%------PRESIONES

load('sin_con_QLh10.mat')%carga datos hs=10 micras
figure(1)
plot(t1s,pres_sinf)
hold on
plot(t,pres_conf)
xlabel('t(s)')
ylabel('Presión (Pa)')
legend('Sin fugas','Con fugas')

axes.ax1=gca; axes.ax1.FontSize=fsize;
figure(2)
%------COMPARATIVA CON DIFERENTES DISTANCIAS ENTRE SLIPPER Y PLATO
%                           (PRESIONES)
plot(t,pres_conf)

load('con_QLh20.mat') %carga datos hs=20micras
hold on
plot(t,pres_conf)  


load('con_QLh25.mat') %carga datos hs=25micras
hold on
plot(t,pres_conf)  

XMIN= 0; XMAX=0.015;YMIN=.2815E8;YMAX=.308E8;
axis([XMIN XMAX YMIN YMAX ])
title('A diferentes distancias entre slipper y plato ')
xlabel('t(s)')
ylabel('Presión de la cámara (Pa)')
legend('h=10 \mum','h=20 \mum','h=25 \mum')
axes.ax2=gca; axes.ax2.FontSize=fsize;

%                               (CAUDALES)

figure(3)
subplot(2,1,1)
load('sin_con_QLh10.mat') 
plot(t,QL1)
XMIN= 0; XMAX=0.015;YMIN=-2.1E-7;YMAX=.5E-7;
axis([XMIN XMAX YMIN YMAX ])
title('QL1')
xlabel('t(s)');
ylabel('Caudal QL1 (m^3 /s)');
axes.ax3=gca; axes.ax3.FontSize=fsize;

subplot(2,1,2)
plot(t,QL2)
title('QL2')
xlabel('t(s)');
ylabel('Caudal QL2 (m^3 /s)');
axes.ax4=gca; axes.ax4.FontSize=fsize;

%------COMPARATIVA CON DIFERENTES DISTANCIAS ENTRE CILINDRO Y PISTÓN
%                           (PRESIONES)

figure (4)
subplot(2,1,2)
plot(t,QL2); hold on
title('QL2 para diferentes distancias entre pistón y cilindro')
xlabel('t(s)');
ylabel('Caudal QL2 (m^3 /s)');


load('con_QLh20')%carga datos hs=20micras
subplot(2,1,2)
plot(t,QL2)

load('con_QLh25') %carga datos hs=25micras
subplot(2,1,2)
plot(t,QL2)
legend('hs=10 \mum','hs=20 \mum','hs=25 \mum');
axes.ax5=gca; axes.ax5.FontSize=fsize;


%COMPARATIVA VARIANDO LA DISTANCIA ENTRE PISTÓN Y CILINDRO 
%(se mantiene hs=10micras)


load('sin_con_QLh10.mat') %carga datos h1=2micras 
figure(4)
subplot(2,1,1);
plot(t,QL1)
axes.ax6=gca; axes.ax6.FontSize=fsize;
figure (5)
plot(t,pres_conf)
axes.ax7=gca; axes.ax7.FontSize=fsize;

fprintf(' ----------- CAUDALES REPRESENTATIVOS')
fprintf(' VARIANDO LA DISTANCIA ENTRE PISTÓN Y CILINDRO -----------')
fprintf('\n\t\t\t\t\t\t\t\t\t((CASO: PISTÓN ASCENDENTE)')
[t,i]=unique(t);
fprintf(' \nSi hs=2 micras:\n - Para t=0, QL1=%.3em^3/s \n',QL1(1));
fprintf(' - Para t=0.0075, QL1=%.3em^3/s',interp1(t,QL1(i),0.0075));

load('con_H3') %carga datos h1=3micras

figure(4)
subplot(2,1,1);
hold on; plot(t,QL1)
figure(5)
hold on; plot(t,pres_conf)

[t,i]=unique(t);
fprintf(' \nSi hs=3 micras:\n - Para t=0, QL1=%.3em^3/s \n',QL1(1));
fprintf(' - Para t=0.0075, QL1=%.3em^3/s',interp1(t,QL1(i),0.0075));

load('con_H4') %carga datos h1=4micras

figure(4); subplot(2,1,1); 
hold on; plot(t,QL1)
figure(5)
hold on
plot(t,pres_conf)

[t,i]=unique(t);
fprintf(' \nSi hs=4 micras:\n - Para t=0, QL1=%.3em^3/s \n',QL1(1));
fprintf(' - Para t=0.0075, QL1=%.3em^3/s',interp1(t,QL1(i),0.0075));

load('con_H5') %carga datos h1=5micras

figure(4); subplot(2,1,1);
hold on; plot(t,QL1)
figure(5)
hold on
plot(t,pres_conf)

[t,i]=unique(t);
fprintf(' \nSi hs=5 micras:\n - Para t=0, QL1=%.3em^3/s \n',QL1(1));
fprintf(' - Para t=0.0075, QL1=%.3em^3/s\n',interp1(t,QL1(i),0.0075));

figure(4)
subplot(2,1,1)
title('QL1 para diferentes distancias entre slipper y plato ');
xlabel('t(s)');
ylabel('Caudal QL1 (m^3 /s)');
legend('h=2 \mum','h=3 \mum','h=4 \mum','h=5 \mum');
figure(5)
title('Presión a diferentes distancias entre pistón y cilindro ')
xlabel('t(s)')
ylabel('Presión de la cámara (Pa)')
legend('h=2 \mum','h=3 \mum','h=4 \mum','h=5 \mum');


%Representación gráfica (CASO:PISTÓN DESCENDENTE)
% Si la presión de aspiración es de 1 bar
load('down_1bar')
figure(6)
plot(t,p_down)
xlabel('t(s)')
ylabel('Presión (Pa)')
title('Presión durante el descenso (Pentrada=1 bar)')
axes.ax8=gca; axes.ax8.FontSize=fsize;


% Si la presión de aspiración es de 7.5 bar 
load('down_7_5bar')
figure(7)
plot(t,p_down)
xlabel('t(s)')
ylabel('Presión (Pa)')
title('Presión durante el descenso (Pentrada=7.5 bar)')
axes.ax9=gca; axes.ax9.FontSize=fsize;


