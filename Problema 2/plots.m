clear all
close all

% Previamente, es necesario tener los datos de las 
% integraciones guardados  si SOLO se desea ver las gráficas,
% teniendo los cálculos hechos. Si se ejecuta el archivo
% principal, este script se ejecuta de forma automática. 


load('Datos_1.mat'); % Carga los datos del cálculo
                     % realizado con el diámetro mayor
figure(1);
subplot(2,1,1);
plot(t,y(:,1));
title('Velocidad del émbolo');
xlabel('Tiempo (s)');
ylabel('Velocidad (m/s)');
subplot(2,1,2);
plot(t,y(:,2));
title('Posición émbolo');
xlabel('Tiempo (s)');
ylabel('Desplazamiento y (m)');
figure(2);
subplot(2,1,1);
plot(t,y(:,3));
title('Aceite');
xlabel('Tiempo (s)');
ylabel('Presión (Pa)');
subplot(2,1,2)
plot(t,y(:,4));
title('Nitrógeno');
xlabel('Tiempo (s)');
ylabel('Presión (Pa)');
fprintf('\nPRESIONES MÁXIMAS OBTENIDAS SEGÚN EL DIÁMETRO DEL ÉMBOLO\n')
fprintf('\t-> Si Dem=%.2f mm, PN2max=%.2f bar, Pacmax=%.2f bar\n',...
    Dem1*1000,max(y(:,4))*1e-5,max(y(:,3))*1e-5)

load('Datos_2.mat'); % Carga los datos del cálculo
                     % realizado con el diámetro intermedio
figure(1);
subplot(2,1,1);
hold on; plot(t,y(:,1)); 

subplot(2,1,2);
hold on; plot(t,y(:,2)); 

figure(2);
subplot(2,1,1);
hold on; plot(t,y(:,3)); 

subplot(2,1,2);
hold on; plot(t,y(:,4)); 

fprintf('\t-> Si Dem=%.2f mm, PN2max=%.2f bar, Pacmax=%.2f bar\n',...
    Dem2*1000,max(y(:,4))*1e-5,max(y(:,3))*1e-5)

fsize=11;
load('Datos_3.mat'); % Carga los datos del cálculo 
                     % realizado con el diámetro menor
figure(1);
subplot(2,1,1);
hold on; plot(t,y(:,1));
lgd.lgd1=legend(['D=' num2str(Dems(1)*10^3) ' mm'],...
    ['D=' num2str(Dems(2)*10^3) ' mm'],['D=' ...
    num2str(Dems(3)*10^3) ' mm']);
axes.ax1=gca; axes.ax1.FontSize=fsize; axes.ax1.XLim=[0 4*pi];
subplot(2,1,2);
hold on; plot(t,y(:,2)); 
lgd.lgd2=legend(['D=' num2str(Dems(1)*10^3) ' mm'],...
    ['D=' num2str(Dems(2)*10^3) ' mm'],['D=' ...
    num2str(Dems(3)*10^3) ' mm']);
axes.ax2=gca; axes.ax2.FontSize=fsize;axes.ax2.XLim=[0 4*pi];

figure(2);
subplot(2,1,1);
hold on; plot(t,y(:,3)); 
lgd.lgd3=legend(['D=' num2str(Dems(1)*10^3) ' mm'],...
    ['D=' num2str(Dems(2)*10^3) ' mm'],['D=' ...
    num2str(Dems(3)*10^3) ' mm']);
axes.ax3=gca; 
axes.ax3.FontSize=fsize;axes.ax3.XLim=[0 4*pi];
subplot(2,1,2);
hold on; plot(t,y(:,4)); 
lgdlgd4=legend(['D=' num2str(Dems(1)*10^3) ' mm'],...
    ['D=' num2str(Dems(2)*10^3) ' mm'],['D=' ...
    num2str(Dems(3)*10^3) ' mm']);
axes.ax4=gca; axes.ax4.FontSize=fsize;axes.ax4.XLim=[0 4*pi]; 

fprintf('\t-> Si Dem=%.2f mm, PN2max=%.2f bar, Pacmax=%.2f bar\n',...
    Dem3*1000,max(y(:,4))*1e-5,max(y(:,3))*1e-5)

load('Datos_2.mat')

figure(3);
subplot(2,2,1); plot(t,y(:,2));
title('Desplazamiento del émbolo');
xlabel('Tiempo (s)'); ylabel('Desplazamiento y (m)');
axes.ax5=gca; axes.ax5.XLim=[0 4*pi];
axes.ax5.FontSize=fsize;
subplot(2,2,2); plot(t,y(:,3));
title('Aceite'); xlabel('Tiempo (s)'); ylabel('Presión (Pa)');
axes.ax6=gca; axes.ax6.XLim=[0 4*pi];
axes.ax6.YLim=[0.5*10^6 2.5*10^6];
axes.ax6.FontSize=fsize; 


subplot(2,2,3); plot(t,y(:,1));
axes.ax7=gca; axes.ax7.XLim=[0 4*pi];axes.ax7.FontSize=fsize;
title('Velocidad del émbolo'); 
xlabel('Tiempo (s)'); 
ylabel('Velocidad (m/s)');
subplot(2,2,4); plot(t,y(:,4));
axes.ax8=gca; axes.ax8.XLim=[0 4*pi]; axes.ax8.FontSize=fsize;
title('Nitrógeno'); xlabel('Tiempo (s)'); ylabel('Presión (Pa)');


