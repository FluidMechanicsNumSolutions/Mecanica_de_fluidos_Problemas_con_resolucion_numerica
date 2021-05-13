%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCRIPT PRINCIPAL DEL PROBLEMA 5 (v1.0.0)

%%%% Los resultados que se obtienen en este programa
%%%% son los asociados al problema 5 del libro:
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

%NOTA: La nomenclatura utilizada es similar a la presentada en el problema.
% Las unidades de cada magnitud física en el Sistema Internacional

% Geometría (CASO A)
rout=0.080;
ro=0.060;
H=0.001;
alpha=70*pi/180;
rext=ro+H;
rin=0.005;

% Otras constantes (CASO A)
a=cos(alpha);
b=log(ro./rin);
c=(H./(2*sin(alpha)))+ro;
D=(rout-ro)/cos(alpha);

% Condiciones de contorno y propiedades físicas (CASO A)
Pin=2e+5;
Pout=1e+5;
mu=857*32*10^(-6); % Viscosidad dinámica del aceite ISO 32

% Cálculo de la presión Po (CASO A)
Po= (Pout+Pin*log((ro+H/(2*sin(alpha))+D*cos(alpha))/(ro+H/...
    (2*sin(alpha))))/(cos(alpha)*log(ro/rin)))/(1+log((ro+H/...
    (2*sin(alpha))+D*cos(alpha))/(ro+H/(2*sin(alpha))))/...
    (cos(alpha)*log(ro/rin)));

% Caudal volumétrico (CASO A)
QtrcilA=-1/6*(Po-Pin)*pi*H^3/(mu*log(ro/rin)); 
fprintf('------> CASO A\n')
fprintf('\t - Caudal circulante: %.2e m^3/s\n',QtrcilA)

%Distribución de presión en la zona cilíndrica (CASO A)
r2=(rin:0.001:ro);
Pc2=Pin+(Po-Pin)*log(r2/rin)*(1/b);
r=[0,r2];
Pcgen=[Pin,Pc2];

%Distribución de presión en la zona troncocónica (CASO A)
x=linspace(0,D,100);
PtrA=Po+(Po-Pin)*log((ro+H/(2*sin(alpha))+x*cos(alpha))/(ro+...
    H/(2*sin(alpha))))/(cos(alpha)*log(ro/rin));

%Fuerza de sustentación sobre superficie cilíndrica (CASO A)
Fpc = Pin*pi*rin^2+Pin*pi*(-rin^2+ro^2)+2*pi*(Po-Pin)*...
    ((1/2)*ro^2*log(ro/rin)+1/4*(rin^2-ro^2))/log(ro/rin);
fprintf(['\t - Fuerza de sustentación sobre'...
    ' superficie cílindrica: %.3f N\n'],Fpc)
fprintf('\t\t-> Debido a la presión: %.3f N\n',Fpc)
fprintf('\t\t-> Debido a los esfuerzos cortantes: %.1f N\n',0)

%Fuerza de sustentación sobre superficie troncocónica (CASO A)
I1=2*Po*pi*(ro*D + D^2/2*cos(alpha))*cos(alpha);
I2p=ro*b*(a*D*log((a*D+c)/c)-a*D+c*log((a*D+c)/c))/a;
I2pp=b*(2*log((a*D+c)/c)*D^2*a^2-D^2*a^2+2*c*D*a-...
    2*c^2*log((a*D+c)/c))/(4*a);
I2=I2p+I2pp;
FptrA=I1+I2;
FctrA=-2*1/2*H*pi*sin(alpha)*(Po-Pin)*(D+H/(2*sin(alpha)*...
    cos(alpha))*log((ro+H*sin(alpha)/2)/(ro+D*cos(alpha)+...
    H*sin(alpha)/2)))/log(ro/rin);
fprintf(['\t - Fuerza de sustentación sobre superficie' ...
    ' troncocónica: %.3f N\n'],FptrA+FctrA)
fprintf('\t\t-> Debido a la presión: %.3f N\n',FptrA)
fprintf('\t\t-> Debido a los esfuerzos cortantes: %.3f N\n',FctrA)

%Fuerza de sustentacion total (CASO A)
FsustotalA=Fpc+FptrA+FctrA;
fprintf(['\t - Fuerza de sustentación TOTAL' ...
    ' (CASO A): %.3f N\n\n'],FsustotalA)

% Geometría (CASO B)
% Se añaden los datos de geometría necesarios para el caso de superficie
% cilíndrica
resf=(ro:0.01:30);
gamma_fin=asin(ro./resf);
gamma_in=atan((rin)./(resf+H));

% Otras constantes (CASO B)
A=cos(alpha);
B2=log(tan(gamma_fin/2)./tan(gamma_in/2));
C=(H./(2*sin(alpha)))+resf;

% Cálculo de la presión P2 
% Nota: Para cada radio de esfera, habrá una presión P2 diferente
P2 = (Pout+Pin.*log((ro+H./(2*sin(alpha))+D*cos(alpha))./...
    (ro+H./(2*sin(alpha))))./(cos(alpha)*log(tan((1/2).*gamma_fin)./...
    tan((1/2).*gamma_in))))./(1+log((ro+H./(2.*sin(alpha))+...
    D.*cos(alpha))./(ro+H./(2*sin(alpha))))./(cos(alpha).*...
    log(tan((1/2).*gamma_fin)./tan((1/2).*gamma_in))));

% Caudal volumétrico (CASO B)
% Nota: Para cada radio de esfera, habrá una caudal circulante distinto
fprintf('------> CASO B\n')
fprintf(['Para cada radio de esfera, el caudal,' ...
    ' la distribución de presión, '])
fprintf(['las fuerzas de sustentación... varían.\n' ...
    ' A continuación, en el Command Window,'])
ind=2;
fprintf([' se presentan los resultados si el radio '...
    'de esfera fuera de %.2f cm \n'],resf(ind)*100)
QtrcilB=-1/6*(P2-Pin)*pi*H^3/(mu*log(tan(gamma_fin/2)/tan(gamma_in/2)));
fprintf('\t - Caudal circulante: %.2e m^3/s\n',QtrcilB(ind))

% Distribución de presión en la zona esférica y troncocónica (CASO B)
% Nota: Para cada radio de esfera, habrá una distribución de presión
% diferente. 
Pesf=zeros(length(P2),2000)+Pin;
PtrB=zeros(length(P2),length(x));
gamma=zeros(length(P2),2000);
for i=1:length(resf)
    gammainicial=gamma_in(i);
    gammafinal=gamma_fin(i);
    gamma1=linspace(0,gammainicial,1000);
    gamma2=linspace(gammainicial,gammafinal,1000);
    gamma(i,1:end)=[gamma1,gamma2];
    B2p=log(tan(gammafinal/2)./tan(gammainicial/2)); 
    Pesf(i,1001:2000)=Pin+((P2(i)-Pin)./B2p).*...
        (log(tan(gamma2./2)/tan(gammainicial/2)));
    PtrB(i,1:length(x))=P2(i)+(P2(i)-Pin).*log((ro+H/(2*sin(alpha))+...
        x*cos(alpha))/(ro+H/(2*sin(alpha))))./(cos(alpha)*...
        log(tan(gamma_fin(i)/2)/tan(gamma_in(i)/2)));
end

% Fuerza de sustentación sobre la superficie esférica 
% debida a la presión (En función del radio de la esfera) 
Fpesf=zeros(1,length(gamma_fin));
indices=[];
for i=1:length(gamma_fin)
    dFpesf= @(x) ((P2(i)-Pin).*log(tan((1./2).*x)./...
        tan((1./2).*gamma_in(i)))./log(tan((1./2).*gamma_fin(i))./...
        tan((1./2).*gamma_in(i)))).*((resf(i)).^2).*pi.*sin(2.*x);
    Fpesf(i)=integral(dFpesf,gamma_in(i),gamma_fin(i))+...
        Pin*pi*rin^2+0.5*Pin*pi*resf(i).^2.*...
        (cos(2*gamma_in(i))-cos(2*gamma_fin(i))); 
    if resf(i)==2 | resf(i)==8 | resf(i)==30
       indices=[indices i];
    end
end

% Fuerza de sustentación sobre la superficie esférica
% debida a los esfuerzos cortantes (En función del radio de la esfera) 
B1=2*pi*(P2-Pin)./B2;
Fcesf=+(B1./2)*H.*(resf.^2).*(1./(resf+H/2)).*...
    (cos(gamma_fin)-cos(gamma_in));
fprintf(['\t - Fuerza de sustentación sobre ' ...
    'superficie esférica: %.3f N\n'],Fpesf(ind)+Fcesf(ind))
fprintf('\t\t-> Debido a la presión: %.4f N\n',Fpesf(ind))
fprintf('\t\t-> Debido a los esfuerzos cortantes: %.3f N\n',Fcesf(ind))

% Fuerza de sustentación sobre la superficie troncocónica (CASO B)
% debida a la presión
I1=2*Po*pi*(ro*D + D^2/2*cos(alpha))*cos(alpha);
I2p=ro*B2.*(A*D*log((A*D+C)./C)-A*D+C.*log((A*D+C)./C))/A;
I2pp=B2.*(2*log((A*D+C)./C)*D^2*A^2-D^2*A^2+...
    2*C*D*A-2*C.^2.*log((A*D+C)./C))/(4*A);
I2=I2p+I2pp;
FptrB=I1+I2;

% Fuerza de sustentación sobre la superficie troncocónica (CASO B)
% debida a los esfuerzos cortantes 
FctrB=-H*pi*sin(alpha)*(P2-Pin).*(D+H/(2*sin(alpha)*cos(alpha))*...
    log((ro+H*sin(alpha)/2)/(ro+D*cos(alpha)+H*sin(alpha)/2)))./...
    log(tan(gamma_fin/2)./tan(gamma_in/2));
fprintf(['\t - Fuerza de sustentación sobre' ...
    ' superficie troncocónica: %.3f N\n'],FptrB(ind)+FctrB(ind))
fprintf('\t\t-> Debido a la presión: %.3f N\n',FptrB(ind))
fprintf('\t\t-> Debido a los esfuerzos cortantes: %.3f N\n',FctrB(ind))

% Fuerza de sustentación total (CASO B)
FsustotalB=Fpesf+Fcesf+FptrB+FctrB;
fprintf(['\t - Fuerza de sustentación TOTAL '...
    '(CASO B): %.3f N\n\n'],FsustotalB(ind))

% Otros resultados que se comentan en el apartado 3 del problema
fprintf(['---- Otros resultados que se comentan' ...
    ' en el apartado 3 del problema \n'])
fprintf('\t==>Fuerza de sustentación sobre superficie esférica debida a ')
fprintf(['la presión\n\t  si los valores de radio ' ...
    'de esfera son elevados:\n'])
for i=indices
    fprintf('\t - Si resf = %.2f m, F = %.3f N\n\n',resf(i),Fpesf(i)); 
end
fprintf('===============SUBA PARA VER MÁS RESULTADOS===============\n')

%Representación gráfica 
fig{1}=figure(1);
fig{1}.Position=[41 309 699 480];
subplot(2,1,1)
plot(r,Pcgen,'LineWidth',1.7)
xlabel('r (m)'); ylabel('Presión (Pa)');
title('Distribución de presión en la zona cilíndrica (CASO A)');
ax=gca;
ax.FontSize=12;

subplot(2,1,2)
plot(x,PtrA,'LineWidth',1.7)
xlabel('x (m)'); ylabel('Presión (Pa)');
title('Distribución de presión en la zona cónica (CASO A)');
ax=gca;
ax.FontSize=12;

fig{2}=figure(2);
fig{2}.Position=[740 306 699 480];
subplot(2,1,1)
plot(gamma(ind,:)*180/pi,Pesf(ind,:),'LineWidth',1.7)
xlabel('Ángulo \gamma (^{o})'); ylabel('Presión (Pa)');
title('Distribución de presión en la zona esférica (CASO B)');
ax=gca;
ax.FontSize=12;

subplot(2,1,2)
plot(x,PtrB(ind,:),'LineWidth',1.7)
xlabel('x (m)'); ylabel('Presión (Pa)');
title('Distribución de presión en la zona cónica (CASO B)');
ax=gca;
ax.FontSize=12;

fig{3}=figure(3);
fig{3}.Position=[211 317 1044 481];
plot(resf,Fpesf,'LineWidth',1.7)
xlabel('Radio de la esfera (m)'); ylabel('Fuerza (N)');
title_name=['Fuerza de sustentación debida a la presión'...
    ' en función del radio de la esfera'];
title(title_name);
ax=gca;
ax.XLim=[0,2];
ax.FontSize=12;

fig{4}=figure(4);
fig{4}.Position=[211 130 1044 481];
plot(resf,Fcesf,'LineWidth',1.7)
xlabel('Radio de la esfera (m)'); ylabel('Fuerza (N)');
title_name=['Fuerza de sustentación debida a los esfuerzos'... 
    ' cortantes en función del radio de la esfera'];
title(title_name);
ax=gca;
ax.XLim=[0,2];
ax.FontSize=12;


