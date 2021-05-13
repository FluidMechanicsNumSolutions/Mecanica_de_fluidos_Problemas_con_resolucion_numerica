function [ dxdt ] = sistodessub( t,x,p )
% En esta funcion se introduce las derivadas de primer orden 
% de las variables que se desean calcular. Para ello es necesario 
% tener los parametros de la geometria, condiciones de contorno 
% y propiedades fisicas. 

% t: tiempo
% x: Variables que se integrarán: 
%    => x(1): Posicion del piston 
%    => x(2): Velocidad del embolo 
%    => x(3): Posicion del embolo
%    => x(4): Presion del aceite en la camara 1
%    => x(5): Presion del aceite en la camara 2
% dxdt: Derivadas de primer orden de cada una de las variables: 
%    => dxdt(1): velocidad del pistón
%    => dxdt(2): aceleración del émbolo 
%    => dxdt(3): velocidad del émbolo
%    => dxdt(4): variación de la presión del aceite en la camara 1 respecto del
%                tiempo
%    => dxdt(5): variación de la presión del aceite en la camara 2 respecto del
%                tiempo
% p: parámetros necesarios definidos en el script principal


Do=p(1);
So=p(2);
Dem=p(3);
Sem=p(4);
St=Sem; 
Dag=p(5);
Sag=p(6);
Ho=p(7);
H1=p(8);
e=p(9); 
Eemb=p(10);
Z=e+H1+Ho;
Ht=p(11);
W=p(12);
A=p(13);
memb=p(14);
g=p(15);
gamma=p(16);
Cd=p(17);
N=p(18);
PN2in=p(19);
Y2in=p(20);
Pacin=p(21);
denin=p(22);
B=p(23);
visc=p(24);
Fc=p(25);
Fs=p(26);
us=p(27);
n=p(28);
Fv=p(29);

y1=x(1); %Posicion del piston 
u2=x(2); %Velocidad del embolo
y2=x(3); %Posicion del embolo
Pac1=x(4); %Presion del aceite en la camara 1
Pac2=x(5); %Presion del aceite en la camara 2

Pac1in=Pacin;
den1=denin*exp((Pac1-Pacin)/B);
den2=denin*exp((Pac2-Pacin)/B);

VN2=St*(Ht-y2-Eemb);
VN2in=St*(Ht-Y2in-Eemb);

reldD=Dag/Do;
m=reldD^2;
relLD=e/Dag;


Ffr=(Fc+(Fs-Fc)*exp((-u2/us)^n)+Fv*u2);


dxdt(1)=A*W*sin(W*t);
u1=dxdt(1);
dxdt(2)=real(-(PN2in*(VN2in/(St*(Ht-Eemb-y2)))^gamma...
    *Sem-Pac2*Sem+memb*g)/memb);
dxdt(3)=u2;

for j=1:6 
    %Bucle para aproximar de forma mas precisa el valor de la 
    %primera derivada teniendo una aproximacion del Cd
    dxdt(4)=real(-(N*Cd*sqrt(2*((Pac1-Pac2)/den1)/(1-(N^2)*(m^2)))*...
        Sag-u1*So)*B/(H1*St+Ho*So-So*y1));   
    
    Qag=(B*So*u1-dxdt(4)*(H1*St+Ho*So-So*y1))/(B*N);
    Re=abs(Dag.*Qag/(visc*Sag));
    Cd=CD2mod2(Re,relLD,m); %La funcion CD2mod2 calcula el Cd
    
    dxdt(5)=real(-(den1*N*Cd*sqrt(2*((Pac1-Pac2)/den1)/(1-(N^2)*(m^2)))*...
        Sag-den2*St*u2)*B/(den2*(St*(-y2+Z))));
 end


dxdt=dxdt(:);
end

