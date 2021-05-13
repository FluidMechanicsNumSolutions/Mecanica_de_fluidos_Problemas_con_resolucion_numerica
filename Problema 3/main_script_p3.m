%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCRIPT PRINCIPAL DEL PROBLEMA 3 (v1.0.0)

%%%% Los resultados que se obtienen en este programa
%%%% son los asociados al problema 3 del libro:
%%%%
%%%% MEC�NICA DE FLUIDOS. PROBLEMAS CON RESOLUCI�N NUM�RICA.
%%%% Autores: �lvaro Ma�as Gonz�lez y Josep M. Bergada Grany� 
%%%%
%%%% Para mayor informaci�n acerca del enunciado del problema 
%%%% y la resoluci�n del mismo, es necesario adquirir el libro completo

%%%% Este script es el principal, el que debe ser ejecutado
%%%% para obtener los resultados definidos en el problema.
%%%% Los datos de entrada pueden ser modificados. 
%%%% Una vez ejecutado este script, se presentan los resultados 
%%%% en Command Window y en las figuras corresponientes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
close all

% Definici�n de la geometr�a
%--> Di�metros y superficies
Do=35*10^(-3); % Di�metro del pist�n (m)
So=pi*(Do^2)/4; %Superficie del pist�n (m^2)
Dem=55*10^(-3); %Di�metro del �mbolo (m)
Sem=pi*(Dem^2)/4; %Superficie del �mbolo (m^2)
St=Sem; %Superficie total (m^2)
Dag=5*10^(-3); %Di�metro del agujero (m)
Sag=pi*(Dag^2)/4; %Superficie del agujero (m^2)

%--> Longitudes definidas en el esquema:
Ho=105*10^(-3);%Longitud Ho (m)
H1=55*10^(-3);%Longitud H1 (m)
e=10*10^(-3); %Longitud agujero (m)
Eemb=5*10^(-3);%Longitud Eemb (m)
Z=e+H1+Ho; %Suma de longitudes (m)
Ht=Z+140*10^(-3); %Longitud total (m)

%--> Relaciones entre di�metros y longitudes
reldD=Dag/Do; %Relaci�n de di�metros
m=reldD^2; %Relaci�n de di�metros al cuadrado
relLD=e/Dag; %Relaci�n entre longitud del agujero y su di�metro

% Definicion de las condiciones de contorno y propiedades f�sicas
freq=20; %Frecuencia del pist�n (Hz)
W=2*pi*freq; %Velocidad angular del pist�n (rad/s)
Amp=2; %Amplitud de la velocidad (m/s)
A=Amp/W; %Amplitud de oscilaci�n del pist�n (m)
memb=0.1;% Masa del embolo (kg)
g=9.81; %Gravedad (m/s^2)
gamma=1.4; % Constante adiab�tica para el nitr�geno
B=1.8*10^9; %Modulo de compresibilidad del aceite (Pa)
visc=42*10^(-6); %Viscosidad cinem�tica (m^2/s)
denac=0.870*1000; %Densidad del aceite (kg/m^3). 
                  %(Es la densidad inicial)
visd=visc*denac; %Viscosidad din�mica (Pa�s)

Cd=0.7; %Coeficiente de descarga.
        %ATENCION!! Este valor es una aproximaci�n inicial.
        %Para cada instante de tiempo se calcular� correctamente.
N=7; %Numero de agujeros
vs=0.1; %Velocidad de Stribeck (m/s)
n=2; %Par�metro experimental
Fs=20; %Fuerza est�tica (N)
Fc=10; %Fuerza de Coulomb (N)
Fv=30; %Coeficiente de fricci�n viscosa (kg/s)

% Definici�n de las condiciones iniciales

PN2in=5*(10^5); %Presi�n inicial del nitr�geno (Pa)
Y2in=Z+60*10^(-3); %Posici�n inicial del embolo (m)
VN2in=St*(Ht-Y2in-Eemb); % Volumen inicial de nitr�geno (m^3)
u1in=0; %Velocidad inicial del pist�n (m/s)
Y1in=40*10^(-3); %Posici�n inicial del pist�n (m)
u2in=0; %Velocidad inicial del embolo (m/s)
Y2in=Z+60*10^(-3);%Posici�n inicial del embolo (m)
Pac1in=(PN2in*Sem+g*memb)/Sem; %Presi�n inicial en la c�mara inferior (Pa)
Pac2in=Pac1in; %Presi�n inicial en la c�mara intermedia (Pa)
denin=denac; %Densidad inicial del aceite en las dos c�maras (kg/m^3)

% Resoluci�n num�rica del sistema
ResolNum; %Se ejecuta el script que resuelve num�ricamente todo el sistema

%--> comprobaci�n
fx=N*Qag-u1*So;

% Representaci�n grafica de la fuerza de amortiguamiento...

represent(t,y1,u1,y2,u2,Pac1,Pac2,PN2,den1,den2,...
    Qag,Rec,Famortcomp,Ff1,Do,false);

[Qmax indice]=max(Qag);
tQmax=t(indice);
[Pmax indice]=max(Pac1);
tPmax=t(indice);
[Famortmax indice]=max(Famortcomp);
tFmax=t(indice);

%--> Para diferentes di�metros del pist�n

y(:,1)=linspace(15,45,4)*10^(-3);
for ii=1:length(y)
    Do=y(ii);
    So=pi*(Do^2)/4;
    reldD=Dag/Do;
    m=reldD^2;
    ResolNum
    z1(ii,:)=Famortcomp;
    z2(ii,:)=Ff1;
    [maxFamort(ii,1) maxFamort(ii,2)]=max(Famortcomp);
    QFmax(ii,1)=Qag(maxFamort(ii,2));
    maxFamort(ii,2)=t(maxFamort(ii,2));
    hold on
end
represent(t,y1,u1,y2,u2,Pac1,Pac2,PN2,den1,den2,...
    Qag,Rec,z1,z2,y,true);

%[[[ Nota en Command Window con alguno de los resultados
% num�ricos comentados en los apartados 3 y 4:
% (Parte del c�digo que permite la visualizaci�n)]]]

fprintf('\n%s\n','APARTADO 3:');
fprintf('\n\t%s %.3e %s \n','- Caudal m�ximo por cada agujero:',...
    Qmax,'m^3/s');
fprintf('\t%s %.4f %s \n','   => Tiempo asociado',tQmax,'s');
fprintf(['\n\t%s %.3e %s \n','- Presi�n m�xima en la c�mara'...
    '1 de aceite:'],Pmax,'Pa');
fprintf('\t%s %.4f %s \n','   => Tiempo asociado',tPmax,'s');
fprintf('\n\t%s %.3f %s \n',...
    '- Fuerza amortiguamiento m�xima en la c�mara 1 de aceite:',...
    Famortmax,'N');
fprintf('\t%s %.4f %s \n','   => Tiempo asociado',tFmax,'s');

fprintf('\n%s\n','APARTADO 4:');
fprintf('\n%s\n',['======================' ...
    'PARA DIFERENTES DI�METROS======================'])

for i=1:length(y)

    fprintf('\n%s %.3f %s \n',...
        '====>>>>> Con un di�metro del pist�n de ',y(i),'metros')
    fprintf('\n\t%s %.2f %s \n',...
        'Fuerza de amortiguamiento m�xima:  ',maxFamort(i,1),' N')
    fprintf('\n\t%s %.4f %s \n',...
        'Tiempo asociado a esta fuerza m�xima:  ',maxFamort(i,2),' s')
    fprintf('\n\t%s %.3e %s \n\t%s \n',...
        'Caudal asociado a esta fuerza m�xima:  ',...
        QFmax(i),' m^3/s','(Por cada agujero)')

end

fprintf('\n\t\t%s\n',['SUBA HACIA ARRIBA PARA VER TODOS' ...
    'LOS RESULTADOS NUM�RICOS COMENTADOS'])
