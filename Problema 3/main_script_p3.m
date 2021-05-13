%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCRIPT PRINCIPAL DEL PROBLEMA 3 (v1.0.0)

%%%% Los resultados que se obtienen en este programa
%%%% son los asociados al problema 3 del libro:
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
clc
close all

% Definición de la geometría
%--> Diámetros y superficies
Do=35*10^(-3); % Diámetro del pistón (m)
So=pi*(Do^2)/4; %Superficie del pistón (m^2)
Dem=55*10^(-3); %Diámetro del émbolo (m)
Sem=pi*(Dem^2)/4; %Superficie del émbolo (m^2)
St=Sem; %Superficie total (m^2)
Dag=5*10^(-3); %Diámetro del agujero (m)
Sag=pi*(Dag^2)/4; %Superficie del agujero (m^2)

%--> Longitudes definidas en el esquema:
Ho=105*10^(-3);%Longitud Ho (m)
H1=55*10^(-3);%Longitud H1 (m)
e=10*10^(-3); %Longitud agujero (m)
Eemb=5*10^(-3);%Longitud Eemb (m)
Z=e+H1+Ho; %Suma de longitudes (m)
Ht=Z+140*10^(-3); %Longitud total (m)

%--> Relaciones entre diámetros y longitudes
reldD=Dag/Do; %Relación de diámetros
m=reldD^2; %Relación de diámetros al cuadrado
relLD=e/Dag; %Relación entre longitud del agujero y su diámetro

% Definicion de las condiciones de contorno y propiedades físicas
freq=20; %Frecuencia del pistón (Hz)
W=2*pi*freq; %Velocidad angular del pistón (rad/s)
Amp=2; %Amplitud de la velocidad (m/s)
A=Amp/W; %Amplitud de oscilación del pistón (m)
memb=0.1;% Masa del embolo (kg)
g=9.81; %Gravedad (m/s^2)
gamma=1.4; % Constante adiabática para el nitrógeno
B=1.8*10^9; %Modulo de compresibilidad del aceite (Pa)
visc=42*10^(-6); %Viscosidad cinemática (m^2/s)
denac=0.870*1000; %Densidad del aceite (kg/m^3). 
                  %(Es la densidad inicial)
visd=visc*denac; %Viscosidad dinámica (Pa·s)

Cd=0.7; %Coeficiente de descarga.
        %ATENCION!! Este valor es una aproximación inicial.
        %Para cada instante de tiempo se calculará correctamente.
N=7; %Numero de agujeros
vs=0.1; %Velocidad de Stribeck (m/s)
n=2; %Parámetro experimental
Fs=20; %Fuerza estática (N)
Fc=10; %Fuerza de Coulomb (N)
Fv=30; %Coeficiente de fricción viscosa (kg/s)

% Definición de las condiciones iniciales

PN2in=5*(10^5); %Presión inicial del nitrógeno (Pa)
Y2in=Z+60*10^(-3); %Posición inicial del embolo (m)
VN2in=St*(Ht-Y2in-Eemb); % Volumen inicial de nitrógeno (m^3)
u1in=0; %Velocidad inicial del pistón (m/s)
Y1in=40*10^(-3); %Posición inicial del pistón (m)
u2in=0; %Velocidad inicial del embolo (m/s)
Y2in=Z+60*10^(-3);%Posición inicial del embolo (m)
Pac1in=(PN2in*Sem+g*memb)/Sem; %Presión inicial en la cámara inferior (Pa)
Pac2in=Pac1in; %Presión inicial en la cámara intermedia (Pa)
denin=denac; %Densidad inicial del aceite en las dos cámaras (kg/m^3)

% Resolución numérica del sistema
ResolNum; %Se ejecuta el script que resuelve numéricamente todo el sistema

%--> comprobación
fx=N*Qag-u1*So;

% Representación grafica de la fuerza de amortiguamiento...

represent(t,y1,u1,y2,u2,Pac1,Pac2,PN2,den1,den2,...
    Qag,Rec,Famortcomp,Ff1,Do,false);

[Qmax indice]=max(Qag);
tQmax=t(indice);
[Pmax indice]=max(Pac1);
tPmax=t(indice);
[Famortmax indice]=max(Famortcomp);
tFmax=t(indice);

%--> Para diferentes diámetros del pistón

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
% numéricos comentados en los apartados 3 y 4:
% (Parte del código que permite la visualización)]]]

fprintf('\n%s\n','APARTADO 3:');
fprintf('\n\t%s %.3e %s \n','- Caudal máximo por cada agujero:',...
    Qmax,'m^3/s');
fprintf('\t%s %.4f %s \n','   => Tiempo asociado',tQmax,'s');
fprintf(['\n\t%s %.3e %s \n','- Presión máxima en la cámara'...
    '1 de aceite:'],Pmax,'Pa');
fprintf('\t%s %.4f %s \n','   => Tiempo asociado',tPmax,'s');
fprintf('\n\t%s %.3f %s \n',...
    '- Fuerza amortiguamiento máxima en la cámara 1 de aceite:',...
    Famortmax,'N');
fprintf('\t%s %.4f %s \n','   => Tiempo asociado',tFmax,'s');

fprintf('\n%s\n','APARTADO 4:');
fprintf('\n%s\n',['======================' ...
    'PARA DIFERENTES DIÁMETROS======================'])

for i=1:length(y)

    fprintf('\n%s %.3f %s \n',...
        '====>>>>> Con un diámetro del pistón de ',y(i),'metros')
    fprintf('\n\t%s %.2f %s \n',...
        'Fuerza de amortiguamiento máxima:  ',maxFamort(i,1),' N')
    fprintf('\n\t%s %.4f %s \n',...
        'Tiempo asociado a esta fuerza máxima:  ',maxFamort(i,2),' s')
    fprintf('\n\t%s %.3e %s \n\t%s \n',...
        'Caudal asociado a esta fuerza máxima:  ',...
        QFmax(i),' m^3/s','(Por cada agujero)')

end

fprintf('\n\t\t%s\n',['SUBA HACIA ARRIBA PARA VER TODOS' ...
    'LOS RESULTADOS NUMÉRICOS COMENTADOS'])
