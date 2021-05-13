
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCRIPT PRINCIPAL DEL PROBLEMA 8 (v1.0.0)

%%%% Los resultados que se obtienen en este programa
%%%% son los asociados al problema 8 del libro:
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

% => La nomenclatura utilizada es la misma que en el problema
% => Las unidades utilizadas están en Sistema Internacional

% Propiedades físicas y condiciones de contorno
Xo=25; 
Yo=45;
rho=1.2;
Gamma=350;
x= linspace(0,90,500); 
y= linspace(0,90,500);

% Distribución de presiones sobre la placa plana
% a lo largo del eje x (ec. 91)
AP_91= distr_P_down(Gamma,rho,x,Xo,Yo);

figure(1); subplot(2,1,1); plot(x,AP_91,'LineWidth',1.2)
xlabel('Posición horizontal (m)'); ylabel('P-P_o (Pa)')
xlim([0 90])
ax=gca; ax.FontSize=10;

% Flujo de vorticidad
Epsilon=diff(AP_91)./diff(x);

subplot(2,1,2); plot(x, [0 Epsilon],'LineWidth',1.2);
xlabel('Posición horizontal (m)'); ylabel('Flujo de vorticidad (Pa/m)')
xlim([0 90])
ax=gca; ax.FontSize=10;

 % Velocidades sobre la placa plana
[Vx_1,Vy_1,V_1]=VxVy(0,Gamma,x,0,Xo,Yo);
 
% Distribución de presiones sobre la placa plana
% a lo largo del eje x (ec. 89)
AP_89=-rho*(V_1.^2)/2; 
 
figure(2); subplot(2,1,1); plot(x,AP_89,'LineWidth',1.2)
xlabel('Posición horizontal (m)'); ylabel('P-P_o (Pa)')
xlim([0 90])
ax=gca; ax.FontSize=10;
 
% Flujo de vorticidad
Epsilon=diff(AP_89)./diff(x);

subplot(2,1,2); 
plot(x,[Epsilon(1) Epsilon],'LineWidth',1.2);
xlabel('Posición horizontal (m)'); ylabel('Flujo de vorticidad (Pa/m)'); 
xlim([0 90])
ax=gca; ax.FontSize=10;

% Evolución a lo largo del tiempo

Gamma_o=Gamma; %Circulación inicial = Circulación definida anteriormente
Te1=125;
tfinal=200;

p=cell(2,4);
Uinf=[0 0.1 0.5 1];
markers=['*';'^';'x';'o'];
t=linspace(0,tfinal,500);
dt=t(2)-t(1);
Uo=zeros(length(t),1); %Inicialización de variable
Vo=zeros(length(t),1); %Inicialización de variable

% ==> Con circulación constante a lo largo del tiempo: 

for i=1:length(Uinf)
    
    % Integración numérica para hallar la posición del vortice en función
    % del tiempo
    [t,VAR]=ode45(@(t,x) UoVo(t,x,Uinf(i),0,Gamma_o,false),t,[Xo, Yo]);
    x_vort=VAR(:,1);
    y_vort=VAR(:,2);
    clear VAR
    Uo=Uinf(i)+(Gamma./(4*pi)).*(x_vort.^2)./...
        (y_vort.*(x_vort.^2 + y_vort.^2));
    Vo=-(Gamma./(4*pi)).*y_vort.^2./(x_vort.*...
        (x_vort.^2 + y_vort.^2));
    
    % => Representación gráfica
    figure(3);
    p{i}=plot(x_vort,y_vort,'LineWidth',1.2); hold on;
    p{i}.Marker=markers(i);
    p{i}.MarkerIndices=fix(logspace(0,log10(length(t)),10));
    xlim([20,120]); ylim([20 50]); 
    xlabel('Posición horizontal (m)'); ylabel('Posición vertical (m)'); 
    ax=gca; ax.FontSize=11;
    figure(4)
    subplot(1,2,1)
    p{i}=plot(t,Uo,'LineWidth',1.2); hold on
    p{i}.Marker=markers(i);
    p{i}.MarkerIndices=fix(logspace(0,log10(length(t)),10));
    xlabel('Tiempo (s)') ; ylabel('Velocidad horizontal (m/s)')
    ax=gca; ax.FontSize=11;
    subplot(1,2,2)
    p{i}=plot(t,Vo,'LineWidth',1.2); hold on
    p{i}.Marker=markers(i);
    p{i}.MarkerIndices=fix(logspace(0,log10(length(t)),10));
    xlabel('Tiempo (s)') ; ylabel('Velocidad vertical (m/s)')
    ax=gca; ax.FontSize=11;
end
legend('U_{\infty}=0 m/s','U_{\infty}=0.1 m/s',...
    'U_{\infty}=0.5 m/s','U_{\infty}=1 m/s');
subplot(1,2,1)
legend('U_{\infty}=0 m/s','U_{\infty}=0.1 m/s',...
    'U_{\infty}=0.5 m/s','U_{\infty}=1 m/s');
figure(3)
legend('U_{\infty}=0 m/s','U_{\infty}=0.1 m/s',...
    'U_{\infty}=0.5 m/s','U_{\infty}=1 m/s');

% ==> Con circulación variable a lo largo del tiempo: 

for i=1:length(Uinf)
    
%%%%% Resolución mediante Euler: 
%%%%% (Nota: Utilizar este código que se introduce a continuación, 
%%%%% o bién, el que no está comentado, es decir, el siguiente  
%%%%% sin comentar (sin %) , es equivalente, los resultados 
%%%%% son muy similares)
%     x_vort=zeros(length(t),1);
%     y_vort=zeros(length(t),1);
%     Uo=zeros(length(t),1);
%     Vo=zeros(length(t),1);
%     x_vort(1)=Xo;
%     y_vort(1)=Yo;
%     
%     for j=1:length(t)-1
%         Gamma=Gamma_o*exp(-t(j)/(Te1));
%         Uo(j)=Uinf(i)+(Gamma/(4*pi))*(x_vort(j)^2)/(y_vort(j)...
%             *(x_vort(j)^2 + y_vort(j)^2));
%         Vo(j)=-(Gamma/(4*pi))*y_vort(j)^2/(x_vort(j)...
%             *(x_vort(j)^2 + y_vort(j)^2));
%         x_vort(j+1)=x_vort(j)+Uo(j)*dt;
%         y_vort(j+1)=y_vort(j)+Vo(j)*dt;
%     end
%     
%%%%% Resolución mediante la función ode45 de MATLAB
%%%%%% Nota: Con ode45 es posible conseguir mayor grado de precisión 
%%%%%% ya que utiliza un método numérico de mayor orden 
    [t,VAR]=ode45(@(t,x) UoVo(t,x,Uinf(i),Te1,Gamma_o,true),t,[Xo, Yo]);
    Gamma=Gamma_o*exp(-t./(Te1));
    x_vort=VAR(:,1);
    y_vort=VAR(:,2);
    clear VAR
    Uo=Uinf(i)+(Gamma./(4*pi)).*(x_vort.^2)./(y_vort.*...
        (x_vort.^2 + y_vort.^2));
    Vo=-(Gamma./(4*pi)).*y_vort.^2./(x_vort.*(x_vort.^2 + y_vort.^2));
%%%%%
    if Uinf(i)==0
        save('dataU0.mat',"x_vort","y_vort","Vo","Uo")
    end
    
    % => Representación gráfica
    figure(5); 
    p{i}=plot(x_vort,y_vort,'LineWidth',1.2); hold on;
    p{i}.Marker=markers(i);
    p{i}.MarkerIndices=fix(logspace(0,log10(length(t)),10));
    xlim([20,100]); ylim([20 50]); 
    xlabel('Posición horizontal (m)'); ylabel('Posición vertical (m)');
    ax=gca; ax.FontSize=11;
    figure(6)
    subplot(1,2,1)
    p{i}=plot(t,Uo,'LineWidth',1.2); hold on
    p{i}.Marker=markers(i);
    p{i}.MarkerIndices=fix(logspace(0,log10(length(t)),10));
    xlabel('Tiempo (s)') ; ylabel('Velocidad horizontal (m/s)')
    ax=gca; ax.FontSize=11;
    subplot(1,2,2)
    p{i}=plot(t,Vo,'LineWidth',1.2); hold on
    p{i}.Marker=markers(i);
    p{i}.MarkerIndices=fix(logspace(0,log10(length(t)),10));
    xlabel('Tiempo (s)') ; ylabel('Velocidad vertical (m/s)')
    ax=gca; ax.FontSize=11;
end
legend('U_{\infty}=0 m/s','U_{\infty}=0.1 m/s',...
    'U_{\infty}=0.5 m/s','U_{\infty}=1 m/s');
subplot(1,2,1)
legend('U_{\infty}=0 m/s','U_{\infty}=0.1 m/s',...
    'U_{\infty}=0.5 m/s','U_{\infty}=1 m/s');
figure(5)
legend('U_{\infty}=0 m/s','U_{\infty}=0.1 m/s',...
    'U_{\infty}=0.5 m/s','U_{\infty}=1 m/s');

clear p

load("dataU0.mat")
[X,Y] = meshgrid(x,y);
ind = fix(linspace(1,length(x),100));

t_sim = 120; % Tiempo final para la representación gráfica
ind_tiempo= find(t>=t_sim,1);
r=zeros(1,7);k=0;

for i=fix(linspace(1,ind_tiempo,5))
    k=k+1;
    % Campo de velocidades en un tiempo = t(i)
    [Vx,Vy,V,r(k)]=VxVy(0,Gamma(i),X,Y,x_vort(i),y_vort(i));
    
    % Distribución de presiones en un tiempo = t(i)
    AP=distr_P_gen(0,Gamma(i),rho,X,Y,x_vort(i),y_vort(i),V);
    
    % => Representación gráfica
    
    f=figure;
    f.Units='normalized';
    f.Position=[0.05 0.56 0.87 0.3];
    ax=gca; ax.FontSize=11;
    
    subplot(1,3,1)
    quiver(X(ind,ind),Y(ind,ind),Vx(ind,ind),Vy(ind,ind),2);
    xlabel('Posición horizontal (m)'); ylabel('Posición vertical (m)');
    title(['Tiempo = ' num2str(round(t(i),0)) ' s'])
    xlim([x_vort(i)-10 x_vort(i)+10]);ylim([y_vort(i)-10 y_vort(i)+10])  
    ax=gca; ax.FontSize=11;
    
    subplot(1,3,2)
    contourf(X,Y,V,linspace(0,25,100),'LineStyle','none');
    xlabel('Posición horizontal (m)'); ylabel('Posición vertical (m)'); 
    title(['Tiempo = ' num2str(round(t(i),0)) ' s'])
    xlim([0 90]); ylim([0 70])
    ax=gca; ax.FontSize=11;
    colormap(jet);
    c=colorbar; c.Label.String='Velocidad (m/s)';
    c.Limits=[0 25];
    shading interp;
    
    subplot(1,3,3)
    contourf(X,Y,AP,linspace(-350,0,100),'LineStyle','none');
    xlabel('Posición horizontal (m)'); ylabel('Posición vertical (m)'); 
    title(['Tiempo = ' num2str(round(t(i),0)) ' s'])
    xlim([0 90]); ylim([0 70])
    ax=gca; ax.FontSize=11;
    colormap(jet);
    c=colorbar; c.Label.String='Presión relativa (Pa)';
    shading interp;
    
end