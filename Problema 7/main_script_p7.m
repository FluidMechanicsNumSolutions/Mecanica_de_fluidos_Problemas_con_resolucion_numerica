%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCRIPT PRINCIPAL DEL PROBLEMA 7 (v1.0.0)

%%%% Los resultados que se obtienen en este programa
%%%% son los asociados al problema 7 del libro:
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

% => La nomenclatura utilizada es la misma que en problema
% => Las unidades utilizadas en las variables están en Sistema Internacional

% Geometría
Rs=[40.1 40.05 40.15]*1e-3;
e=[50 30 10 5 1]*1e-6;
R=0.040;
hp=Rs-R;
psi=hp/R;

%Condición de contorno y propiedades físicas
Pin=0*(10.^6);
w=1000*pi;
mu=0.02; 
Pcav=0;

% Otras variables 
n=1000; % Número de puntos por cada caso
fignum=0; %Variable para establecer una referencia en figuras
          % (utilizada en la representación gráfica)
Pentrada=zeros(length(Rs),length(e)); % Inicialización de variable
Pentrada_v=zeros(length(Rs),n); % Inicialización de variable
angin=zeros(length(Rs),length(e)); % Inicialización de variable
ang_pin=zeros(length(Rs),n); % Inicialización de variable
e_vect=zeros(length(Rs),n); % Inicialización de variable
Fsu2=zeros(length(Rs),n); % Inicialización de variable
L2=zeros(length(Rs),n); % Inicialización de variable
fig=cell(1,15);

%Resolución de los casos
for j=1:length(Rs)
    fprintf('>>>>>>>>>> CASO %i\n',j)
    %Resolución de integrales 
    anginicial=-2*pi/9; % Ángulo inicial
    angfinal=16*pi/9; % Ángulo final
    anggen=linspace(anginicial,angfinal,n); % Ángulo genérico
    w_vect=linspace(0,1000*pi,n);
    P=zeros(length(anggen),length(e));
    Fsu=zeros(length(w_vect),length(e));
    i2p=zeros(1,length(anggen));
    i3p=zeros(1,length(anggen));
    norepresent=0;
    
    for k=1:length(e)
        epsi=e(k)/hp(j);
        if abs(epsi-1)<0.001
            norepresent=k;
            Pentrada(j,k)=NaN;
        else
            di2dx=@(x) (1-epsi.*cos(x)).^(-2); 
            di3dx=@(x) (1-epsi.*cos(x)).^(-3);
            di4dx= @(x) ((1-epsi*cos(x)).^(-2)).*cos(x);
            di5dx= @(x) ((1-epsi*cos(x)).^(-3)).*cos(x);
            for i= 1:length(anggen)
                i2p(i)=integral(di2dx,anginicial,anggen(i));
                i3p(i)=integral(di3dx,anginicial,anggen(i));
            end
            i2=integral(di2dx,anginicial,angfinal);
            i3=integral(di3dx,anginicial,angfinal);
            i4=integral(di4dx,anginicial,angfinal);
            i5=integral(di5dx,anginicial,angfinal);

            fprintf('--> Si e=%.3f mm:\n',e(k)*1e+3)
            fprintf('\ti2=%.4f\ti3=%.4f\t',i2,i3)

            % Distribución de presión en función del ángulo
            P(1:length(i2p),k)=Pin+(6*mu*w/psi(j)^2)*(i2p-(i2./i3)*(i3p));

            % Presión de entrada obtenida a partir de la gráfica
            Pentrada(j,k)=Pcav+abs(min(P(1:length(i2p),k))); 
                
            fprintf('\tP%i%imin=%.3e Pa\t',k,j,min(P(1:length(i2p),k))) 

            % Fuerza de sustentación en función de la velocidad angular
            Fsu(1:length(w_vect),k)=-(6*mu*R*w_vect/psi(j)^2)*(i4-i2*i5/i3);
            
            % Longitud del cojinete en función de la velocidad angular
            L1=18000./(2.*Fsu);
            fprintf('\tL%i%i=%.3f cm\n',k,j,interp1(w_vect,L1(1:length(w_vect),k)*100,w))
        end

    end
    ex=e;
    
    % Fuerza de sustentación en función de la excentricidad
    e_vect(j,1:n)=linspace(0.00005*hp(j),hp(j)-0.00005*hp(j),n);
    epsi=e_vect(j,:)/hp(j);
    di2dx=@(x) (1-epsi.*cos(x)).^(-2); 
    di3dx=@(x) (1-epsi.*cos(x)).^(-3); 
    di4dx=@(x) ((1-epsi.*cos(x)).^(-2)).*cos(x);
    di5dx=@(x) ((1-epsi.*cos(x)).^(-3)).*cos(x);
    ii2=integral(di2dx,anginicial,angfinal,'ArrayValued',true);
    ii3=integral(di3dx,anginicial,angfinal,'ArrayValued',true);
    ii4=integral(di4dx,anginicial,angfinal,'ArrayValued',true);
    ii5=integral(di5dx,anginicial,angfinal,'ArrayValued',true);

    Fsu2(j,1:length(e_vect))=-(6*mu*R*w/psi(j)^2)*(ii4-ii2.*ii5./ii3);
    
    % Longitud del cojinete en función de la excentricidad
    L2(j,1:length(e_vect))=18000./(2.*Fsu2(j,1:length(e_vect)));

    % Presión de entrada en función de la excentricidad
    
    ang_pin(j,1:length(e_vect(j,:)))=acos((1-ii2./ii3)./epsi);
    d2Pdang2=6*mu*w*(-2*epsi.*sin(ang_pin)./(1-epsi.*cos(ang_pin)).^3+...
        3*ii2.*epsi.*sin(ang_pin)./(ii3.*(1-epsi.*cos(ang_pin)).^4))/psi(j).^2;
    if isequal(d2Pdang2<0, zeros(length(Rs),n)) 
        di2dx=@(x) (1-epsi.*cos(x)).^(-2); 
        di3dx=@(x) (1-epsi.*cos(x)).^(-3);
        for i=1:length(ang_pin)
             i2p(i)=integral(@(x) (1-epsi(i)*cos(x)).^(-2),anginicial,ang_pin(j,i));
             i3p(i)=integral(@(x) (1-epsi(i)*cos(x)).^(-3),anginicial,ang_pin(j,i));
        end
        Pentrada_v(j,1:length(e_vect(j,:)))=Pcav+abs((6*mu*w/psi(j)^2)...
            *(i2p-(ii2./ii3).*(i3p)));
    
    else
       disp('No se cumple el criterio de mínimo relativo en al menos uno de los ángulos') 
    end
    % Número de Sommerfeld en función de la excentricidad
    cte=(psi(j)^2)/(mu*w*R);
    S=Fsu2(j,1:length(e_vect))*((psi(j)^2)/(mu*w*R));
    S2=6*(-ii4+ii2.*ii5./ii3);

    %Representación gráfica

    cont=0;
    fignum=fignum+2;
    fig{fignum-1}=figure(fignum-1);
    fig{fignum-1}.Position=[9+450*(j-1) 400 452 361];
    fig{fignum}=figure(fignum);
    fig{fignum}.Position=[9+450*(j-1) 375 452 361];
    
    for i=1:k
        if norepresent==0 
            leg{i,1}=num2str(ex(i)*1000);
            figure(fignum-1)
            p{1}(i)=plot(anggen,P(:,i),'LineWidth',1.25); hold on;
            figure(fignum)
            p{2}(i)=plot(w_vect,Fsu(:,i),'LineWidth',1.25); hold on;
        else
            if i==norepresent
            else
                cont=cont+1;
                leg{cont,1}=num2str(ex(i)*1000);
                figure(fignum-1)
                p{1}(i)=plot(anggen,P(:,i),'LineWidth',1.25); hold on;
                figure(fignum)
                p{2}(i)=plot(w_vect,Fsu(:,i),'LineWidth',1.25); hold on;
            end
        end
    end
    
    p{1}(1).Marker='o'; p{2}(1).Marker='o';
    p{1}(2).Marker='*'; p{2}(2).Marker='*';  
    p{1}(3).Marker='^'; p{2}(3).Marker='^';  
    p{1}(4).Marker='s'; p{2}(4).Marker='s';  
    p{1}(5).Marker='d'; p{2}(5).Marker='d';  
    
    np=5;
    
    p{1}(1).MarkerIndices=(1:n/np:n); 
    p{2}(1).MarkerIndices=(1:n/np:n); 
    p{1}(2).MarkerIndices=(fix(0.5*n/np):n/np:n);
    p{2}(2).MarkerIndices=(fix(0.5*n/np):n/np:n);  
    p{1}(3).MarkerIndices=(1:n/np:n); 
    p{2}(3).MarkerIndices=(1:n/np:n);  
    p{1}(4).MarkerIndices=(fix(0.5*n/np):n/np:n); 
    p{2}(4).MarkerIndices=(fix(0.5*n/np):n/np:n);  
    p{1}(5).MarkerIndices=(1:n/np:n); 
    p{2}(5).MarkerIndices=(1:n/np:n); 
    
    figure(fignum-1)
    xlim([-pi/4 7*pi/4]); xticks(-pi/4:pi/4:7*pi/4)
    xticklabels({'-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4',...
        '3\pi/2','7\pi/4'})
    xlabel('\psi (rad)'); ylabel('Presión (Pa)')
    legend(leg)
    ax=gca; ax.FontSize=11;

    figure(fignum);
    
    xlim([0 1000*pi]); xticks(0:125*pi:1000*pi)
    xticklabels({'0\pi','125\pi','250\pi','375\pi','500\pi','625\pi',...
        '750\pi','875\pi','1000\pi'})
    xlabel('\omega (rad/s)'); ylabel('Fuerza de sustentación (N/m)')
    legend(leg)
    ax=gca; ax.FontSize=11;
    

    fignum=fignum+1;
    fig{fignum}=figure(fignum);
    fig{fignum}.Position=[9+450*(j-1) 350 452 361];
    
    plot(e_vect(j,:)*1000,Fsu2(j,1:length(e_vect)),'LineWidth',1.5)
    xlabel('e (mm)'); ylabel('Fuerza de sustentación (N/m)')
    ax=gca; ax.FontSize=11; ax.XLim=[0 max(e_vect(j,:))*1000+0.01];
    ax.YLim=[-0.5*1e+8 max(Fsu2(j,1:length(e_vect)))]; 

    fignum=fignum+1;
    fig{fignum}=figure(fignum);
    fig{fignum}.Position=[9+450*(j-1) 325 452 361];

    plot(e_vect(j,:)*1000,S,'LineWidth',1.5)
    xlabel('e (mm)'); ylabel('Número de Sommerfeld')
    ax=gca; ax.FontSize=11; ax.XLim=[0 max(e_vect(j,:))*1000+0.01]; 
    ax.YLim=[-100 max(S)]; 
    clear leg
    
end

fignum=fignum+1;
fig{fignum}=figure(fignum);
fig{fignum}.Position=[50 150 450 360];
    
for i=1:length(Rs)
    p{2}(i)=plot(e_vect(i,:)*1000,Pentrada_v(i,:),'LineWidth',1.5); hold on
end
p{2}(1).LineStyle='--'; 
p{2}(2).LineStyle='-.';
p{2}(3).LineStyle=':';

xlabel('e (mm)'); ylabel('Presión entrada (Pa)')
xlim([0 max(hp)*1000]); ylim([0 max(Pentrada,[],'all')]);
legend('Caso 1','Caso 2', 'Caso 3')
ax=gca; ax.FontSize=11;

fignum=fignum+1;
fig{fignum}=figure(fignum);
fig{fignum}.Position=[500 150 450 360];

for i=1:length(Rs)
    p{3}(i)=plot(e_vect(i,:)*1000,L2(i,:)*100,'LineWidth',1.5); hold on
end
p{3}(1).LineStyle='--'; 
p{3}(2).LineStyle='-.';
p{3}(3).LineStyle=':';

xlabel('e (mm)'); ylabel('Longitud del cojinete (cm)')
xlim([-0.005 min(hp)*1000]); ylim([-5 max(L2(:,1))]);
legend('Caso 1','Caso 2', 'Caso 3')
ax=gca; ax.FontSize=11;

fignum=fignum+1;
fig{fignum}=figure(fignum);
fig{fignum}.Position=[950 150 450 360];

for i=1:length(Rs)
    p{4}(i)=plot(e_vect(i,:)*1000,Fsu2(i,:),'LineWidth',1.5); hold on
end
p{4}(1).LineStyle='--'; 
p{4}(2).LineStyle='-.';
p{4}(3).LineStyle=':';

xlabel('e (mm)'); ylabel('Longitud del cojinete (cm)')
xlim([-0.005 min(hp)*1000]); ylim([-5 max(L2,[],'all')]);
xlabel('e (mm)'); ylabel('Fuerza de sustentación (N/m)')
ax=gca; ax.FontSize=11; ax.XLim=[0 max(hp)*1000+0.01]; 
ax.YLim=[-0.25*1e+8 min(Fsu2(:,end))]; 
legend('Caso 1','Caso 2', 'Caso 3')
