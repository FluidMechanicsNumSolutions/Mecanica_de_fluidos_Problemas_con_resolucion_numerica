function [ dxdt ] = sistodesbaj( t,x,p)
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
Pac1in=p(21);
Pac2in=p(22);
denin1=p(23);
denin2=p(24);
B=p(25);
visc=p(26);
Fc=p(27);
Fs=p(28);
us=p(29);
n=p(30);
Fv=p(31);

y1=x(1); %Posicion del piston 
u2=x(2); %Velocidad del embolo
y2=x(3); %Posicion del embolo
Pac1=x(4); %Presion del aceite en la camara 1
Pac2=x(5); %Presion del aceite en la camara 2

reldD=Dag/Do;
m=reldD^2;
relLD=e/Dag;

den1=denin1*exp((Pac1-Pac1in)/B);
den2=denin2*exp((Pac2-Pac2in)/B);

VN2=St*(Ht-y2-Eemb);
VN2in=St*(Ht-Y2in-Eemb);

Ffr=(-Fc-(Fs-Fc)*exp((-u2/us)^n)+Fv*u2);

dxdt(1)=A*W*sin(W*t+pi);
u1=dxdt(1);
dxdt(2)=real(-(PN2in*(VN2in/(St*(Ht-Eemb-y2)))^gamma...
    *Sem-Pac2*Sem+memb*g)/memb);
dxdt(3)=u2;

for j=1:6

    dxdt(4)=real(+(den2*N*Cd*sqrt(2*((Pac2-Pac1)/den2)...
        /(1-(N^2)*(m^2)))*Sag+den1*u1*So)*B/(den1*(H1*St+Ho*So-So*y1)));

    Qag=-(B*So*u1-dxdt(4)*(H1*St+Ho*So-So*y1))/(B*N);
    Re=abs(Dag.*Qag/(visc*Sag));
    Cd=CD2mod2(Re,relLD,m);

    dxdt(5)=real(+(N*Cd*sqrt(2*((Pac2-Pac1)/den2)...
        /(1-(N^2)*(m^2)))*Sag+St*u2)*B/(St*(-y2+Z)));

end

dxdt=dxdt(:);
end

