%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCRIPT PRINCIPAL DEL PROBLEMA 9 (v1.0.0)

%%%% Los resultados que se obtienen en este programa
%%%% son los asociados al problema 9 del libro:
%%%%
%%%% MECÁNICA DE FLUIDOS. PROBLEMAS CON RESOLUCIÓN NUMÉRICA.
%%%% Autores: Álvaro Mañas González y Josep M. Bergada Granyó 
%%%%
%%%% Para mayor información acerca del enunciado del problema 
%%%% y la resolución del mismo, es necesario adquirir el libro completo

%%%% Este script es el principal, el que debe ser ejecutado
%%%% para obtener los resultados definidos en el problema.
%%%% Los datos de entrada pueden ser modificados. 

%%%% Para cada uno de los casos analizados en el problema, se debe
%%%% introducir los datos de entrada asociados a cada caso estudiado.
%%%% En Command Window, se presentan los resultados según las variables 
%%%% de entrada introducidas

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
clc
clear all

global L12 L34 D1 D3 Sup1 Sup3 rE12 rE34 R Po1 To1 gamma viscosdin

%%%% Introducción de datos
% Geometría
L12=6; %Longitud de 1 a 2 (m)
L34=10; %Longitud de 3 a 4 (m)
D1=0.05; % Diámetro 1 (m)
D2=D1;
D3=0.1; % Diámetro 3 (m) 
Sup1=pi*(D1^2)/4; % Superficie 1 (m^2)
Sup3=pi*(D3^2)/4; % Superficie 2 (m^2)
E12=0.008*10^(-3);% Rugosidad absoluta entre 1 y 2 (m)
E34=0.008*10^(-3);% Rugosidad absoluto entre 3 y 4 (m)
rE12=E12/D1; % Rugosidad relativa entre 1 y 2
rE34=E34/D3; % Rugosidad absoluto entre 3 y 4

% Propiedades físicas y condiciones de contorno
R=287; %Constante de los gases (J/kg*K)
Po1=4*10^5; %presion de estancamiento en 1 (Pa)
To1=370; %temperatura de estancamiento en 1 (K)
gamma=1.4; %Índice de politropía
viscosdin=2.181*(10^(-5)); %viscosidad dinamica del aire (Pa*s)
P4real=1*10^5;

To2=To1;
To3=To2;
To4=To3;

if Po1>P4real
%%%% 1- Determinar si el flujo está bloqueado en el punto 2
% BAJO LA HIPÓTESIS: NÚMERO DE MACH EN 2 M2=1.

M2=1;
% Cálculo del número de Mach en el punto 1:

%Para el cálculo del número de Mach en 1, se realiza un bucle con el
%objetivo de hallar el número de Reynolds adecuado y consecuentemente, el
%factor de fricción que más se ajuste. De esta manera, la hipótesis que se
%plantea en la resolución del enunciado, donde se supone un número de 
%Reynolds elevado, se comprueba directamente en el inicio del programa. 

M1=0.5; % Suposición inicial para el numero de Mach 1
Re1=massflow(M1,Po1,To1,Sup1,R,gamma)*D1/(Sup1*viscosdin);
i=0;
f1=frictionfact(Re1,rE12);
M1new=FannoMachNum(f1*L12/D1,gamma);
if length(M1new)>1
    M1new=M1new(1);
end

while abs(M1new-M1)>0.00001
    M1=M1new;
    Re1=massflow(M1,Po1,To1,Sup1,R,gamma)*D1/(Sup1*viscosdin);
    f1=frictionfact(Re1,rE12);
    M1new=FannoMachNum(f1*L12/D1,gamma);
    if length(M1new)>1
        M1new=M1new(1);
    end
    i=i+1;
end
M1=M1new;

% Calculado M1, es posible hallar la presión y temperatura 
% en el punto 1 y 2. 

T1=To1/relToT(M1,gamma);
P1=Po1/relPoP(M1,gamma);
P2=P1/Fannorel2P(M1,gamma);
Po2=Po1/Fannorel2Po(M1,gamma);
T2=T1/Fannorel2T(M1,gamma);
To2=To1;

% Cálculo del número de Mach en 3:

relA2A3=Sup3/Sup1; 
M3=IsentMachNum(relA2A3,M2,gamma);
%Se obtienen dos soluciones para M3. Con esto se obtendran dos posibles
%resultados para el punto 3 y 4. Se descartara uno de ellos

Po3=Po2; 
To3=To2; 

% Cálculo del número de Mach en 4:
% También se calcula la temperatura y presión en el punto 3
M4=[];
for i=1:length(M3)
    M=M3(i);
    fLD3=FannoIndvar(M,gamma);
    Re3=massflow(M3(i),Po3,To3,Sup3,R,gamma)*D3/(Sup3*viscosdin);
    f34=frictionfact(Re3,rE34);
    varindep=(fLD3-(f34*L34/D3));
    M4c=FannoMachNum(varindep,gamma);
    % ¡Atención! En el caso que no sea posible obtener un numero de Mach en
    % 4 para uno de los números de Mach calculados en 3, el resultado de la
    % anterior función será NaN. 
    if length(M4c)>1
        for k=1:2
            if M<1 & M4c(k)<1
                M4=[M4 M4c(k)];
            end
            if M>1 & M4c(k)>1
                M4=[M4 M4c(k)];
            end
        end
    else
        M4=[M4 M4c]; 
    end

end

P1=P1+zeros(1,length(M4));
P2=P2+zeros(1,length(M4));

T1=T1+zeros(1,length(M4));
T2=T2+zeros(1,length(M4));
To4=To3;
for i=1:length(M4)
    T3(i)=To3/relToT(M3(i),gamma);
    P3(i)=Po3/relPoP(M3(i),gamma);
    T4(i)=To4/relToT(M4(i),gamma);
end

M1=M1+zeros(1,length(M4));
M2=M2+zeros(1,length(M4));
M3=M3';
% Evaluación de relaciones de presión entre 1 y 4 (Po1/P4):

r14=zeros(1,length(M4));
for i=1:length(M4)
    if M4(i)~=0
        r14(i)=(Po1/Po2)*(Po2/Po3)*Fannorel2Po(M3(i),gamma)*...
            (1/Fannorel2Po(M4(i),gamma))*relPoP(M4(i),gamma);
    end
end
P4=Po1./r14;
if (M4(1)-1)<0.01
    fluj4bloqueado=true; %
else
    fluj4bloqueado=false;
end
fluj24supsonic=false; %Asignación inicial.  
                      %Si el flujo entre 2 y 4 fuera supersónico, 
                      %se determinará más adelante. 
casomax=1+length(M4)*0.1;
CASO=[1.1:0.1:casomax];

ondchoq=zeros(1,length(CASO)); %Para indicar si se tiene en cuenta  
%la onda de choque. En estos casos calculados no se ha tenido en 
%cuenta por ello se rellena con ceros. Posteriormente, 
%los demás casos si se tendrá en cuenta la onda de choque. 
longsinond=length(M4); %Longitud del vector una vez realizado 
                       %el cálculo sin onda de choque
indOK=0; %Índice utilizado guardar saber la posición del caso posible
            % - si indOK=0: inicialización de la variable o caso NO
            %               encontrado
            % - si indOK diferente de 0: posición del caso posible

clc

%%%% Se añaden estos resultados en la tabla de resultados mostrada en el
%%%% 'Command Window'.  
comentario=''; %Variable para añadir cualquier comentario adicional

fprintf('%s',['======================================'...
    '================================'])
fprintf('%s',' HIPÓTESIS 1 M2=1 ')
fprintf('%s \n',['===================================='...
    '=================================='])
fprintf('%s \n','===========>>>>>> Sin onda de choque')
fprintf('%7s %14s %7s %7s %7s %7s ','CASO',...
    'Onda choque','M1','M2','M3','M4')
fprintf('%8s %7s %9s %8s %9s %7s %7s %7s ',...
    'P1','P2','P3','P4','T1','T2','T3','T4')
fprintf('%8s %10s %8s \n','r14','r14(real)','Posible')

for i=1:longsinond
     
    if ondchoq(i)==0
        onda{i}='NO';
    else onda{i}='SI';
    end
    if abs(r14(i)-Po1/P4real)>0.01 | r14(i)==0 | isnan(r14(i))
        if M3(i)>1 & M4(i)>1
            pos{i}='--'; %significa que más adelante se comprovará si 
            % realmente es supersónico entre 2 y 4 o no lo es. 
        else
            if fluj4bloqueado
                pos{i}='SI';
                indok=i;
            else
                pos{i}='NO';
            end
        end
    else
        pos{i}='SI';
        indok=i;
    end
    fprintf('%7.1f %10s %13.4f %5.0f %9.4f %7.4f ',...
        CASO(i),onda{i},M1(i),M2(i),M3(i),M4(i))
    if isnan(P4(i))
       fprintf('%9.2e %7.2e %7.2e %5.2e ',P1(i),P2(i),P3(i),P4(i)) 
       fprintf('%11.2f %7.2f %7.2f %7.2f ',T1(i),T2(i),T3(i),T4(i))
       fprintf('%7.3f %8.2f %9s \n',r14(i),Po1/P4real,pos{i})
    else
       fprintf('%9.2e %7.2e %7.2e %7.2e ',P1(i),P2(i),P3(i),P4(i))
       fprintf('%8.2f %7.2f %7.2f %7.2f ',T1(i),T2(i),T3(i),T4(i))
       fprintf('%7.3f %8.2f %9s \n',r14(i),Po1/P4real,pos{i})
    end
    
end

for i=1:length(r14)
    if r14(i)<=Po1/P4real & ne(isnan(r14(i)),1)
        fprintf('\n%s \n',['El flujo está bloqueado en el ' ...
            'punto 2 ==> HIPÓTESIS 1 CORRECTA'])
        HIPOTESIS1=true; 
        break
    else
        if r14(i)>Po1/P4real & ne(isnan(r14(i)),1) & M3(i)<1 & M4(i)<1
            fprintf('\n%s','El flujo NO está bloqueado en el punto 2')
            fprintf('%s \n',[' ==> HIPÓTESIS 1 NO ES '...
                'CORRECTA ==> El flujo es SUBSÓNICO'])
            HIPOTESIS1=false;
            break
        end
    end
     if isnan(r14(1))
         HIPOTESIS1=false;
     end
end

if HIPOTESIS1==true %HIPOTESIS 1: El flujo está bloqueado en 2

if indOK==0
%%%% 2- Si la onda de choque se produjera en el punto 3

for i=1:length(M3)
    if M3(i)>1 & M3(i)~=M3(i-1)
        
% Cálculo del número de Mach aguas abajo de la onda de choque:

        M3bc=MNumNormShock(M3(i),gamma);
        M3b=M3bc;
% Cálculo del número de Mach en en punto 4:
%También se calculan las presiones y temperaturas

        fLD3b=FannoIndvar(M3b,gamma);
        To3b=To3;
        T3b=T3(i)*r2TNormShock(M3(i),gamma);
        Po3=Po2;
        Po3b=Po3*r2PoNormShock(M3(i),gamma);
        P3b=Po3b/relPoP(M3b,gamma);
        Re3b=massflow(M3b,Po3b,To3b,Sup3,R,gamma)*D3/(Sup3*viscosdin);
        f34b=frictionfact(Re3b,rE34);
        varindep=(fLD3b-(f34b*L34/D3));
        M4c=FannoMachNum(varindep,gamma);
        % De esta última función se obtienen dos soluciones, M4c>1 y
        % M4c<21. Después de una onda de choque, el flujo es subsónico y si
        % es flujo es de Fanno, como máximo se podrá bloquear pero nunca
        % superar la unidad. La solución supersónica se descarta en la
        % tabla de resultados añadida en 'Command Window'
           
        for k=1:length(M4c)
            if M4c(k)<1
                M3=[M3 M3(i)];
                P3=[P3 P3(i)];
                T3=[T3 T3(i)];
                M4=[M4 M4c(k)];
                T4=[T4 To4/relToT(M4c(k),gamma)];
            end
            if isnan(M4c)
                M3=[M3 M3(i)];
                P3=[P3 P3(i)];
                T3=[T3 T3(i)];
                M4=[M4 M4c(k)];
                T4=[T4 To4/relToT(M4c(k),gamma)];
            end
        end
    end
end

% Evaluación de relaciones de presión entre 1 y 4 (Po1/P4):

longconond3=length(M4); %longitud del vector M4 una vez incluido 
                        %el cálculo con onda de choque en 3
for i=longsinond+1:length(M3)
    if M4(i)~=0
        for k=1:length(M3)
            if M3(k)>1
            M3sup=M3(k);
            end
        end
        r14(i)=(Po1/Po2)*(Po2/Po3)*(1/r2PoNormShock(M3sup,gamma))*...
            Fannorel2Po(M3b,gamma)*(1/Fannorel2Po(M4(i),gamma))...
            *relPoP(M4(i),gamma);
    end
end

P1=P1(1)+zeros(1,length(M4));
P2=P2(1)+zeros(1,length(M4));
P4=[Po1./r14];

T1=T1(1)+zeros(1,length(M4));
T2=T2(1)+zeros(1,length(M4));

M1=M1(1)+zeros(1,length(M4));
M2=M2(1)+zeros(1,length(M4));

casomax=2+(longconond3-longsinond)*0.1;
CASO=[CASO 2.1:0.1:casomax];
ondchoq=[ondchoq 1+zeros(1,length(M4c))];

%%%% Se añaden estos resultados en la tabla de resultados mostrada en el
%%%% 'Command Window'. 

fprintf('\n%s \n','===========>>>>>> Con onda de choque en el punto 3')
fprintf('%7s %14s %7s %7s %7s %9s %9s ','CASO',...
    'Onda choque','M1','M2','M3','M3b','M4')
fprintf('%8s %8s %9s %9s %7s ','P1','P2','P3','P3bis','P4')
fprintf('%8s %7s %7s %8s %6s ','T1','T2','T3','T3bis','T4')
fprintf('%8s %11s %8s \n','r14','r14(real)','Posible')

for i=longsinond+1:longconond3

    if ondchoq(i)==0
        onda{i}='NO';
    else onda{i}='SI';
    end
    if abs(r14(i)-Po1/P4real)>0.01 | r14(i)==0 | isnan(r14(i))
        pos{i}='NO';
        if abs(r14(i)-Po1/P4real)<0.1
            pos{i}='NO*';
        end
    else
        pos{i}='SI';
        indOK=length(pos);
    end

    fprintf('%7.1f %10s %13.4f %5.0f %9.4f %9.4f %9.4f ',CASO(i),...
        onda{i},M1(i),M2(i),M3(i),M3b,M4(i))
    fprintf('%9.2e %7.2e %7.2e %7.2e %7.2e ',P1(i),P2(i),P3(i),P3b,P4(i))
    fprintf('%8.2f %7.2f %7.2f %7.2f %7.2f ',T1(i),T2(i),T3(i),T3b,T4(i))
    fprintf('%7.3f %7.3f %8s \n',r14(i),Po1/P4real,pos{i})
end

if ((r14(end)>Po1/P4real & abs(r14(end)-Po1/P4real)>0.001 )| isnan(M4c))
    fprintf('\n %s\n',...
        'Es posible que haya una onda de choque entre el punto 2 y 3')
    fprintf('\n %s\n','(Calculando.....)')

%%%% 3- Si la onda de choque se produce en la tobera (entre 2 y 3)


Dt=D2+(D3-D2)/2; % Diámetro inicial estimado
Dant=D2;% Extremo inferior inicial==> Diámetro inferior
Dpos=D3;% Extremo superior inicial==> Diámetro superior
okerror=false;
k=1;

while okerror==false
    
% Cálculo del número de Mach aguas arriba y aguas abajo de la onda de
% choque en el punto estimado:
% También se calcula temperaturas y presiones
    A2A23=((Dt)^2)/((D2)^2);
    M23=IsentMachNum(A2A23,1,gamma); 
    % M23 solo puede ser supersónico. Con la siguiente condición se
    % selecciona el valor supersonico de la solución calculada. 
    if M23(1)>1
        M23=M23(1);
    elseif M23(2)>1
        M23=M23(2); 
    end
    Po23=Po2;
    P23=Po23/relPoP(M23,gamma);
    Po23b=Po23*r2PoNormShock(M23,gamma);
    To23=To2;
    T23=To23/relToT(M23,gamma);
    To23b=To23;
    Po3p=Po23b;
    To3p=To23b;
    M23b=MNumNormShock(M23,gamma); %Solo se obtiene una solución de esta ec.
    P23b=Po23b/relPoP(M23b,gamma);
    T23b=To23b/relToT(M23b,gamma);
    
% Cálculo del número de Mach en el punto 2 y en el punto 3:
% También se calcula temperaturas y presiones
    A23A3=((D3)^2)/((Dt)^2);
    M3c=IsentMachNum(A23A3,M23b,gamma); %Aparecen dos soluciones
    %M23b es subsónico después de la onda de choque.Por tanto,
    %si aumenta la sección, el número de Mach disminuye, no puede aumentar.
    %Se selecciona cual es el valor de M3 siguiendo este criterio.  
    if M3c(1)<M23b
        M3c=M3c(1);
    elseif M3c(2)<M23b
        M3c=M3c(2); 
    else
        Dcritico=Dt;
    end
    T3p=To3p/relToT(M3c,gamma);
    P3p=Po3p/relPoP(M3c,gamma);
    Re3p=massflow(M3b,Po3p,To3p,Sup3,R,gamma)*D3/(Sup3*viscosdin);
    f34p=frictionfact(Re3p,rE34);
    fLD3p=FannoIndvar(M3c,gamma);
    varindep=(fLD3p-(f34p*L34/D3));
    M4c=FannoMachNum(varindep,gamma); %Pueden aparecer dos soluciones. 
    % Se selecciona el valor que esté por debajo de 1 (subsónico)
    if length(M4c)>1
        if M4c(1)<1
            M4c=M4c(1);
        elseif M4c(2)<1
            M4c=M4c(2);
        else Dt=Dant; break;
        end
    end
% Evaluación de la relación de presión teórica entre el punto 1 y el punto 4:

    r14p=(Po1/Po2)*(Po2/Po23)*(Po23/Po23b)*(Po23b/Po3p)...
        *Fannorel2Po(M3c,gamma)*(1/Fannorel2Po(M4c,gamma))...
        *relPoP(M4c,gamma);
    if abs(r14p-Po1/P4real)/(Po1/P4real)<0.001
        okerror=true;
        P4p=Po1/r14p;
        Po4=P4p*relPoP(M4c,gamma);
    else
        if r14p<Po1/P4real
            Dant=Dt;
            Dt=Dt+(Dpos-Dt)/2;
        else
            Dpos=Dt;
            Dt=Dant+(Dt-Dant)/2;
        end 
        if abs(M4c-1)<0.01 | k>50 
            %El flujo puede quedar bloqueado en el punto 4. 
            okerror=true;
            fluj4bloqueado=true;
        else
            fluj4bloqueado=false;
        end
    end
k=k+1;
end
r14=[r14 r14p];
P3=[P3 P3p];
T3=[T3 T3p];
P4=[Po1./r14];
M3=[M3 M3c];
M4=[M4 M4c];
T4=[T4 To4/relToT(M4c,gamma)];

P1=P1(1)+zeros(1,length(M4));
P2=P2(1)+zeros(1,length(M4));

T1=T1(1)+zeros(1,length(M4));
T2=T2(1)+zeros(1,length(M4));

M1=M1(1)+zeros(1,length(M4));
M2=M2(1)+zeros(1,length(M4));

CASO=[CASO 3];
ondchoq=[ondchoq 1];

fprintf('\n%s','===========>>>>>> Con onda de choque entre 2 y 3');
fprintf(' %10s%0.5f %s \n','con diÁmetro de secciÓn calculado igual a '...
    ,Dt,'metros');
fprintf('%7s %14s %7s %7s  %7s %9s %8s %9s ','CASO',...
    'Onda choque','M1','M2','M23','M23b','M3','M4');
fprintf('%8s %7s %9s %10s %7s %7s ','P1','P2','P23','P23b','P3','P4');
fprintf('%8s %7s %8s %7s %6s %7s ','T1','T2','T23','T23b','T3','T4');
fprintf('%9s %11s %8s \n','r14','r14(real)','Posible');

    if ondchoq(end)==0
        onda{end+1}='NO';
    else onda{end+1}='SI';
    end
    if (abs(r14(end)-Po1/P4real)>0.1 | r14(end)==0 | isnan(r14(end)))...
            & ne(fluj4bloqueado,true)
        pos{end+1}='NO';
    else
        pos{end+1}='SI';
        indOK=length(pos);
        
    end
    fprintf('%7.1f %10s %13.4f %5.0f %9.4f %9.4f %9.4f %9.4f ',...
        CASO(end),onda{end},M1(end),M2(end),M23,M23b,M3(end),M4(end));
    fprintf('%9.2e %7.2e %7.2e %7.2e %7.2e %7.2e ',P1(end), P2(end),...
        P23,P23b,P3(end),P4(end));
    fprintf('%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f ',T1(end), T2(end),...
        T23,T23b,T3(end),T4(end));
    fprintf('%8.3f %8.3f %8s \n',r14(end),Po1/P4real,pos{end});
    

    
%%%% 4- Si la onda de choque se produjera entre los puntos 3 y 4

else
    if pos{end}=='NO'
        fprintf('\n %s \n',...
            'Es posible que haya onda de choque entre el punto 3 y 4')
        fprintf('\n %s\n','(Calculando.....)')

% Determinar en qué posición del conducto comprendido entre 3 y 4 se
% encuentra la onda de choque:

    
for i=1:length(M3)
    if M3(i)>1 & M3(i)~=M3(i-1)
        L334=L34/2; %Suposicion inicial
        L334inf=0; %Respecto el punto 3, la posición del extremo inferior
                   %inicial es nulo.  
        L334sup=L34;%Respecto el punto 3, el extremo superior inicial 
                    %es igual a la longitud del tramo 3-4.
        r14c=-999*1e+10; %Valor imposible (solo para iniciar el bucle)
        k=1;
        M4c=-1; %Valor imposible (solo para iniciar el bucle)
        fluj4bloqueado=false;
        fluj24supsonic=false;
 while (abs(r14c-Po1/P4real)>0.001 | isnan(r14c)) & (M4c<=1 | ...
         isnan(M4c))& k<100 & ne(fluj4bloqueado,true) & ...
         ne(fluj24supsonic,true)
     
     
% Cálculo del número de Mach aguas arriba y aguas abajo de la onda de
% choque en el punto estimado:

        fL3=FannoIndvar(M3(i),gamma);
        %f34=0.0115;
        varindep=(fL3-(f34*L334/D3));
        M34c=FannoMachNum(varindep,gamma);
    if isnan(M34c)
        L334sup=L334;
        L334=(L334+L334inf)/2;
    else
        if length(M34c)>1
            M34=M34c(find(M34c>1,1));
        else
            M34=M34c;
        end
                
        M34bc=MNumNormShock(M34,gamma);
        M34b=M34bc;
        
% Cálculo del número de Mach en el punto 4:

        fL34b=FannoIndvar(M34b,gamma);
        To34=To3;
        T34=To34/relToT(M34,gamma);
        T34b=T34*r2TNormShock(M34,gamma);
        To34b=To34;
        Po34=Po3*Fannorel2Po(M34,gamma)*(1/Fannorel2Po(M3(i),gamma));
        P34=Po34/relPoP(M34,gamma);
        Po34b=Po34*r2PoNormShock(M34,gamma);
        P34b=Po34b/relPoP(M34b,gamma);
        Re34b=massflow(M34b,Po34b,To34b,Sup3,R,gamma)*D3/(Sup3*viscosdin);
        f34b=frictionfact(Re34b,rE34);
        varindep=(fL34b-(f34b*(L34-L334)/D3));
        M4c=FannoMachNum(varindep,gamma);
        if isnan(M4c)
            L334sup=L334;
            L334=(L334+L334inf)/2;
        else
        if length(M4c)>1
            if M4c(1)<1
                M4c=M4c(1);
            elseif M4c(2)<1
                M4c=M4c(2);
            end
        end
        
% Evaluación de la relación de presión teórica 
% entre el punto 1 y el punto 4:

         r14c=(Po1/Po2)*(Po2/Po3)*(Po3/Po34)*(Po34/Po34b)*...
            Fannorel2Po(M34b,gamma)*(1/Fannorel2Po(M4c,gamma))*...
            relPoP(M4c,gamma);
        if abs(r14c-Po1/P4real)>0.001 & M4c<=1 
            if r14c>Po1/P4real
                L334sup=L334;
                L334=(L334+L334inf)/2;
                
            else
                L334inf=L334;
                L334=(L334+L334sup)/2;
                
            end
        end
        if abs(1-M4c)<0.005 & r14c<Po1/P4real & ...
                abs(L334sup-L334inf)*100/L334inf<0.001
            fluj4bloqueado=true;
        end
        if abs(L334-L34)*100/L334inf<0.001 
            %El flujo puede ser supersónico desde el punto 2 hasta el 4
            fluj24supsonic=true;
            
        end
        end
        
    end
 k=k+1;
 end
        M4=[M4 M4c];
        T4=[T4 To4/relToT(M4c,gamma)];
        r14=[r14 r14c];
        T3=[T3 To4/relToT(M3(i),gamma)];
        M3=[M3 M3(i)];
        P3=[P3 P3(i)];
        
    end
end
P1=P1(1)+zeros(1,length(M4));
P2=P2(1)+zeros(1,length(M4));
P4=[Po1./r14];

T1=T1(1)+zeros(1,length(M4));
T2=T2(1)+zeros(1,length(M4));

M1=M1(1)+zeros(1,length(M4));
M2=M2(1)+zeros(1,length(M4));

CASO=[CASO 3];
ondchoq=[ondchoq 1];

fprintf('\n%s','===========>>>>>> Con onda de choque entre 3 y 4')
fprintf([' %10s%0.3f %s \n','a una distancia '...
    'respecto del punto 3 igual a '],L334,'metros')
fprintf('%7s %14s %7s %7s  %7s %9s %8s %9s ','CASO',...
    'Onda choque','M1','M2','M3','M34','M34b','M4');
fprintf('%8s %7s %9s %10s %7s %7s ','P1','P2','P3','P34','P34b','P4');
fprintf('%8s %7s %8s %7s %6s %7s ','T1','T2','T3','T34','T34b','T4');
fprintf('%9s %11s %8s \n','r14','r14(real)','Posible');
    if ondchoq(end)==0
        onda{end+1}='NO';
    else
        onda{end+1}='SI';
    end
    if (abs(r14(end)-Po1/P4real)>0.01 | r14(end)==0 | isnan(r14(end)))...
            & ne(fluj4bloqueado,true) & ne(fluj24supsonic,false)
        pos{end+1}='NO';
    else
        pos{end+1}='SI';
        indOK=length(pos);
    end
   
    if fluj24supsonic & abs(r14(2)-Po1/P4real)>0.01
        
       if r14(2)>Po1/P4real
           comentario=[comentario ...
               'Immediatamente despues al punto 4, se produce '];
           comentario=[comentario 'una compresión brusca en el fluido'];
       else
           if r14(2)<Po1/P4real
                comentario=[comentario ...
                    'Immediatamente despues al punto 4, se produce '];
                comentario=[comentario ...
                    'una expansión brusca en el fluido'];
           end
        end
    end

    fprintf('%7.1f %10s %13.4f %5.0f %9.4f %9.4f %9.4f %9.4f ',...
        CASO(end),onda{end},M1(end),M2(end),M3(end),M34,M34b,M4(end));
    fprintf('%9.2e %7.2e %7.2e %7.2e %7.2e %7.2e ',P1(end), P2(end),...
        P3(end),P34,P34b,P4(end));
    fprintf('%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f ',T1(end), T2(end),...
        T3(end),T34,T34b,T4(end));
    fprintf('%8.3f %8.3f %8s \n',r14(end),Po1/P4real,pos{end});
    

if fluj24supsonic==true
    CASO=[CASO 4];
    ondchoq=[ondchoq 0];
    pos{2}='SI';
    indOK=2;
    fprintf('\n%s',...
        '===========>>>>>> No hay onda de choque entre los puntos 3 y 4. ')
    fprintf('%s \n','El flujo es supersónico entre los puntos 2 y 4.')
    fprintf('%7s %14s %7s %7s %7s %9s ',...
        'CASO','Onda choque','M1','M2','M3','M4')
    fprintf('%8s %7s %9s %8s %9s %7s %7s %7s ',...
        'P1','P2','P3','P4','T1','T2','T3','T4')
    fprintf('%8s %10s %8s \n','r14','r14(real)','Posible')
    fprintf('%7.1f %10s %13.4f %7.4f %7.4f %7.4f ',...
        CASO(2),onda{2},M1(2),M2(2),M3(2),M4(2))
    fprintf('%9.2e %7.2e %7.2e %5.2e ',P1(2),P2(2),P3(2),P4(2)) 
    fprintf('%11.2f %7.2f %7.2f %7.2f ',T1(2),T2(2),T3(2),T4(2))
    fprintf('%7.3f %8.2f %9s \n',r14(2),Po1/P4real,pos{2})
   
end
        
        
    end

end
end
else
    
%%%% 5- Si el flujo es subsónico


M2it=0.5; %Suposición inicial para iniciar el bucle
M2inf=0;
M2sup=1;
[r14c,M1s,M3c,M4c,P1s,P2s,P3s]=Subsonico(M2it);
i=1;

while (abs(r14c(1)-Po1/P4real)>0.001 & i<100 | ...
        ne(isnan(r14c(1)),false)) & abs((M2sup-M2inf)/M2inf)>0.0001
    if isnan(r14c(1))
            M2sup=M2it;
            M2it=(M2inf+M2it)/2;
    else
       if r14c(1)>Po1/P4real
            M2sup=M2it;
            M2it=(M2inf+M2it)/2;
       else
           if r14c(1)<Po1/P4real
               M2inf=M2it;
               M2it=(M2sup+M2it)/2;
           end
       end
    end
    [r14c,M1s,M3c,M4c,P1s,P2s,P3s,P4s,T1s,T2s,T3s,T4s]=Subsonico(M2it);
    
    i=i+1;
end
if abs(1-M4c(1))<0.01
    fluj4bloqueado=true;
else
    fluj4bloqueado=false
end
r14=[r14 r14c(1)];
M1=[M1 M1s];
M2=[M2 M2it];
M3=[M3 M3c(1)];
M4=[M4 M4c(1)];

T1=[T1 T1s];
T2=[T2 T2s];
T3=[T3 T3s];
T4=[T4 T4s];
    
P1=[P1 P1s];
P2=[P2 P2s];
P3=[P3 P3s(1)];
P4=[Po1./r14];

casomax=2+(length(M4c(1)))*0.1;
CASO=[CASO 2.1:0.1:casomax];
ondchoq=[ondchoq zeros(1,length(CASO))];
   
fprintf('\n%s',['=================================' ...
    '====================================='])
fprintf('%s',' HIPÓTESIS 2 M2<1 ')
fprintf('%s \n',['=================================' ...
    '====================================='])

fprintf('%s \n','===========>>>>>> Subsónico en todo el conducto')
fprintf('%7s %14s %7s %7s %7s %7s ',...
    'CASO','Onda choque','M1','M2','M3','M4')
fprintf('%8s %7s %9s %8s %9s %7s %7s %7s ',...
    'P1','P2','P3','P4','T1','T2','T3','T4')
fprintf('%8s %10s %8s \n','r14','r14(real)','Posible')

for i=longsinond+1:length(CASO)
    if ondchoq(i)==0
        onda{i}='NO';
    else onda{i}='SI';
    end
    if abs(r14(i)-Po1/P4real)>0.01 | r14(i)==0 | isnan(r14(i))
        if fluj4bloqueado
            pos{i}='SI';
            indOK=i;
        else
            pos{i}='NO';
        end
    else
        pos{i}='SI';
        indOK=i;
    end
    fprintf('%7.1f %10s %13.4f %7.4f %7.4f %7.4f ',...
        CASO(i),onda{i},M1(i),M2(i),M3(i),M4(i))
    fprintf('%9.2e %7.2e %7.2e %5.2e ',P1(i),P2(i),P3(i),P4(i))
    fprintf('%8.2f %7.2f %7.2f %7.2f ',T1(i),T2(i),T3(i),T4(i))
    fprintf('%7.3f %8.2f %9s \n',r14(i),Po1/P4real,pos{i})
    
end
% 
% 
end

if fluj4bloqueado  & abs(r14(end)-Po1/P4real)>0.01
           comentario='El flujo está bloqueado en 4. ';
       if r14(end)>Po1/P4real
           comentario=[comentario ...
               'Immediatamente despues al punto 4, se produce '];
           comentario=[comentario ...
               'una compresión brusca en el fluido.'];
       else
           if r14(end)<Po1/P4real
                comentario=[comentario ...
                    'Immediatamente después al punto 4, se produce '];
                comentario=[comentario ...
                    'una expansión brusca en el fluido.'];
           end
       end
end
fprintf('\nEl caso posible es el %s. \n%s \n',...
    num2str(CASO(indOK)),comentario)
fprintf('El caudal másico circulante es: %.2f kg/s \n',...
    massflow(M1(indOK),Po1,To1,Sup1,R,gamma))

else
    if Po1==P4real
         fprintf(['No hay diferencia de presión entre' ...
             ' los dipósitos, no puede haber flujo entre ellos \n'])
    end
    if P4real>Po1
         fprintf(['ERROR: Este programa ha sido definido'...
             ' para flujos desde el punto 1 al 4 , NO al revés \n'])
         fprintf(['\t   Por tanto, la presión en 1 debe'...
             ' ser mayor que en 4 para que el flujo vaya'...
             ' en el sentido programado \n'])
    end
end

