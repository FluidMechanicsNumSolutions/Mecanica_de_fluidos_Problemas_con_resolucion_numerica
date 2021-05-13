%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCRIPT PRINCIPAL DEL PROBLEMA 10 (v1.0.0)

%%%% Los resultados que se obtienen en este programa
%%%% son los asociados al problema 10 del libro:
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


%%%% Introducción de datos
% => La nomenclatura utilizada se corresponde a la del problema 
% => Las unidades utilizadas están en Sistema Internacional

% Geometría

X=[5,1.5,1,1,1.5]; % Longitudes
alpha=[8,12,25,30]*pi/180; % Ángulo de inclinación del sólido 

% Propiedades físicas y condiciones de contorno

R=287; % Constante del aire
gamma=1.4; % Constante politrópica

Po=1197; % Presión del aire a una altura de 30000 m
To=226.5; % Temperatura del aire a una altura de 30000 m

Mi=7; % Número de Mach aguas arriba

%%%%% APARTADOS 1,2,3

fprintf(['\n=========================== APARTADOS 1,2,3'...
    ' ===========================\n'])

[M,P,T,F,Fx,Fy,drag,lift,power] = x_43_fun(X,Mi,alpha,gamma,Po,To,R);

%%%%% NOTA: En Command Window se mostrarán los resultados calculados hasta
%%%%%       ahora

nombre_variables=["i","M_i","P_i(Pa)","T_i(K)"];
tabla_condiciones_fluido=table((1:5)',...
    M,P,T,'VariableNames',nombre_variables);
fprintf('\n -> Condiciones de velocidad, presión y temperatura:\n\n')
disp(tabla_condiciones_fluido)

nombre_variables_2=["i","F_i (N)","Fx_i (N)","Fy_i (N)"];
tabla_fuerzas=table((1:5)',F,Fx,Fy,'VariableNames',nombre_variables_2);
fprintf('\n -> Fuerzas sobre las superficies:\n\n')
disp(tabla_fuerzas)

fprintf('\t >>> Fuerza de arrastre:     %.0f N\n',drag)
fprintf('\t >>> Fuerza de sustentación: %.0f N\n',lift)
fprintf('\t >>> Potencia necesaria:     %.0f kW\n',power/1000)

fprintf(['====================================' ...
    '===================================\n']);
disp(' Programa ejecutándose ... ')
%%%%% APARTADO 4

n=50; % Número de puntos 

% >>>> Variando la velocidad de la aeronave
Mi_vect=linspace(2,10,n);

%Inicialización de variables
drag_vect=zeros(n,1); 
lift_vect=drag_vect; 
power_vect=drag_vect; 
M_matrix=zeros(n,5);
P_matrix=zeros(n,5);

for i=1:n
    
    [M,P,T,F,Fx,Fy,drag,lift,power] = ...
        x_43_fun(X,Mi_vect(i),alpha,gamma,Po,To,R);
    drag_vect(i)=drag;
    lift_vect(i)=lift;
    power_vect(i)=power;
    
    for j=1:5
        M_matrix(i,j)=M(j);
        P_matrix(i,j)=P(j);
    end
end

markers=["-o","-^","-d","-+","-"]; 
markerind=fix(linspace(1,n,20));
f=cell(1,2);

for i=1:5
    
    f{1}=figure(1);
    plot(Mi_vect,M_matrix(:,i),markers(i),...
        'MarkerSize',4.5,'MarkerIndices',markerind,'LineWidth',1.2);
    hold on
    
    f{2}=figure(2);
    plot(Mi_vect,P_matrix(:,i),markers(i),...
        'MarkerSize',4.5,'MarkerIndices',markerind,'LineWidth',1.2); 
    hold on
    
end

figure(1)

xlabel('M_{\infty}'); ylabel('Número de Mach'); 
legend('M_1','M_2','M_3','M_4','M_5'); ax=gca;
ax.FontSize=11; 


figure(2)

xlabel('M_{\infty}'); ylabel('Presión (Pa)'); 
legend('P_1','P_2','P_3','P_4','P_5'); ax=gca;
ax.FontSize=11; ylim([-0.5,5]*1e+4);

f{1}.Units='Normalized'; f{1}.Position=[0.15 0.4 0.35 0.4];
f{2}.Units='Normalized'; f{2}.Position=[0.5 0.4 0.35 0.4];


figure
subplot(2,1,1)
plot(Mi_vect,drag_vect,'-.','LineWidth',1.2)
hold on

plot(Mi_vect,lift_vect,'-','LineWidth',1.2)  
hold on

xlabel('M_{\infty}'); ylabel('Fuerza (N)'); 
legend('F. Arrastre','F. Sustentación')
ylim([0, 4*1e+4]); 

subplot(2,1,2)
plot(Mi_vect,power_vect/1000,'LineWidth',1.2)
xlabel('M_{\infty}'); ylabel('Potencia (kW)');

pause(0.1); 
% >>>> Variando el ángulo alpha 1
% (El resto de ángulos se mantienen)

alpha_new=alpha;
alpha_1_min = atan((((gamma+1)./2.*Mi.^2./(Mi.^2.*...
    (1.001/7).^2-1)-1).*tan(asin(1.001/7)))^-1);
alpha1_vect=linspace(alpha_1_min*180/pi,20,n)*pi/180; 

for i=1:n
    
    alpha_new(1)=alpha1_vect(i);    
    [M,P,T,F,Fx,Fy,drag,lift,power] = ...
        x_43_fun(X,Mi,alpha_new,gamma,Po,To,R);
    drag_vect(i)=drag;
    lift_vect(i)=lift;
    power_vect(i)=power;
    
    for j=1:5
        M_matrix(i,j)=M(j);
        P_matrix(i,j)=P(j);
    end
    
end

markers=["-o","-^","-d","-+","-"]; 
markerind=fix(linspace(1,n,20));

for i=1:5

    f{1}=figure(4);
    plot(alpha1_vect*180/pi,M_matrix(:,i),markers(i),...
        'MarkerSize',4.5,'MarkerIndices',markerind,'LineWidth',1.2);
    hold on
   
    f{2}=figure(5);
    plot(alpha1_vect*180/pi,P_matrix(:,i),markers(i),...
        'MarkerSize',4.5,'MarkerIndices',markerind,'LineWidth',1.2); 
    hold on
    
end

figure(4)

xlabel('\alpha_{1} (º)'); ylabel('Número de Mach'); 
legend('M_1','M_2','M_3','M_4','M_5'); ax=gca;
ax.FontSize=11; xlim([0 alpha1_vect(end)*180/pi])

figure(5)

xlabel('\alpha_{1} (º)'); ylabel('Presión (Pa)'); 
legend('P_1','P_2','P_3','P_4','P_5'); ax=gca;
ax.FontSize=11; ylim([-0.5,5]*1e+4);
xlim([0 alpha1_vect(end)*180/pi])

f{1}.Units='Normalized'; f{1}.Position=[0.15 0.4 0.35 0.4];
f{2}.Units='Normalized'; f{2}.Position=[0.5 0.4 0.35 0.4];

figure

subplot(2,1,1)
plot(alpha1_vect*180/pi,drag_vect,'-.','LineWidth',1.2)
hold on

plot(alpha1_vect*180/pi,lift_vect,'-','LineWidth',1.2)  
hold on
xlabel('\alpha_{1} (º)'); ylabel('Fuerza (N)'); 
legend('F. Arrastre','F. Sustentación')
ylim([0, 4*1e+4]); xlim([0 alpha1_vect(end)*180/pi])

subplot(2,1,2)
plot(alpha1_vect*180/pi,power_vect/1000,'LineWidth',1.2)
xlabel('\alpha_{1} (º)'); ylabel('Potencia (kW)');
xlim([0 alpha1_vect(end)*180/pi])

pause(0.1)

% >>>> Variando los ángulos 2 y 3, manteniendo el grosor Y
% (El resto de ángulos se mantienen)

% Cálculos previos
Y=(X(2)*tan(alpha(2))+X(3)*tan(alpha(3)));

alpha_2_max=atan(Y/(X(2)+X(3))); % Ángulo alpha 2 máximo posible. 
epsilon2_min = asin(1.001/7);
alpha_2_min = atan((((gamma+1)./2.*Mi.^2./(Mi.^2.*...
    sin(epsilon2_min).^2-1)-1).*tan(epsilon2_min))^-1);

alpha2_vect=linspace(alpha_2_min,alpha_2_max-1e-5,n); 
        % Obsérvese que al ángulo máximo se le ha restado un valor
        % suficientemente pequeño para que posteriomente el caso límite sea
        % calculable. 
alpha_new=alpha;

% Cálculo de los distintos casos para diferentes angulos alpha 2 y 3
% (El resto de ángulos se mantienen)

for i=1:n

    alpha_new(2)=alpha2_vect(i);
    alpha_new(3)=atan((Y-tan(alpha_new(2))*X(2))/X(3));
    [M,P,T,F,Fx,Fy,drag,lift,power] = ...
        x_43_fun(X,Mi,alpha_new,gamma,Po,To,R);
    drag_vect(i)=drag;
    lift_vect(i)=lift;
    power_vect(i)=power;
  
    for j=1:5
        M_matrix(i,j)=M(j);
        P_matrix(i,j)=P(j);
    end
   
end

markers=["-o","-^","-d","-+","-"]; 
markerind=fix(linspace(1,n,20));

for i=1:5
    
    f{1}=figure(7);
    plot(alpha2_vect*180/pi,M_matrix(:,i),markers(i),...
        'MarkerSize',4.5,'MarkerIndices',markerind,'LineWidth',1.2);
    hold on
    
    f{2}=figure(8);
    plot(alpha2_vect*180/pi,P_matrix(:,i),markers(i),...
        'MarkerSize',4.5,'MarkerIndices',markerind,'LineWidth',1.2); 
    hold on
    
end

figure(7)

xlabel('\alpha_{2} (º)'); ylabel('Número de Mach'); 
legend('M_1','M_2','M_3','M_4','M_5'); ax=gca;
ax.FontSize=11; xlim([0 alpha_2_max*180/pi])

figure(8)

xlabel('\alpha_{2} (º)'); ylabel('Presión (Pa)'); 
legend('P_1','P_2','P_3','P_4','P_5'); ax=gca;
ax.FontSize=11; ylim([-0.5,5]*1e+4);
xlim([0 alpha_2_max*180/pi])

f{1}.Units='Normalized'; f{1}.Position=[0.15 0.4 0.35 0.4];
f{2}.Units='Normalized'; f{2}.Position=[0.5 0.4 0.35 0.4];

figure
subplot(2,1,1)

plot(alpha2_vect*180/pi,drag_vect,'-.','LineWidth',1.2)
hold on

plot(alpha2_vect*180/pi,lift_vect,'-','LineWidth',1.2)
hold on

xlabel('\alpha_{2} (º)'); ylabel('Fuerza (N)');
legend('F. Arrastre','F. Sustentación')
ylim([0, 4*1e+4]); xlim([0 alpha_2_max*180/pi])

subplot(2,1,2)
plot(alpha2_vect*180/pi,power_vect/1000,'LineWidth',1.2)
xlabel('\alpha_{2} (º)'); ylabel('Potencia (kW)');
xlim([0 alpha_2_max*180/pi])

disp(' Fin del programa ')
disp([' SUBA HACIA ARRIBA EN COMMAND WINDOW PARA VER LOS '...
    'RESULTADOS DE LOS APARTADOS 1,2 y 3 '])
