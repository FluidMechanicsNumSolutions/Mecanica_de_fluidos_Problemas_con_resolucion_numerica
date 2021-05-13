%
%
% Este script incorpora toda la resolucion numerica a partir de los datos
% dados en el script principal
%
%
%Script ( 'ResolNum.m' ) debe ser guardado con este nombre junto con los
%demas archivos que componen todo el codigo
%
tfinal=1/freq; 


% Parametros necesarios para la funcion ode23s 
param(1)=Do;
param(2)=So;
param(3)=Dem;
param(4)=Sem;
param(5)=Dag;
param(6)=Sag;
param(7)=Ho;
param(8)=H1;
param(9)=e; 
param(10)=Eemb;
param(11)=Ht;
param(12)=W;
param(13)=A;
param(14)=memb;
param(15)=g;
param(16)=gamma;
param(17)=Cd;
param(18)=N;
param(19)=PN2in;
param(20)=Y2in;
param(21)=Pac1in;
param(22)=denin;
param(23)=B;
param(24)=visc;
param(25)=Fc;
param(26)=Fs;
param(27)=vs;
param(28)=n;
param(29)=Fv;

CI(1)=Y1in;
CI(2)=u2in;
CI(3)=Y2in;
CI(4)=Pac1in;
CI(5)=Pac2in;

%%%===>     SUBIDA DEL PISTON
tfin=(1/2)*(tfinal); 

%Calculo de las variables Y1,Y2, V1,V2, P1, P2 en la subida
%NOTA: Pac1 y Pac2 es lo mismo que P1,P2(nomenclatura en el analisis)
% Si no se definen los errores relativos y absolutos, la funcion ode23s
% trabaja con los errores por defecto. Error rel: 1e-3. Error abs:1e-6.

%options=odeset('RelTol',1e-6,'AbsTol',1e-9);
[ts,ys]=ode23s(@(ts,x) sistodessub(ts,x,param),[0 tfin], CI);
y1s=ys(:,1);
u1s=A*W*sin(W*ts);
u2s=ys(:,2);
y2s=ys(:,3);
Pac1s=ys(:,4);
Pac2s=ys(:,5);

%%%%Recalculo del resto de magnitudes para la subida

VN2s=(St*(Ht-Eemb-y2s));
PN2s= PN2in.*((VN2in./VN2s).^(gamma));
den1s=denin.*exp((Pac1s-Pac1in)./B);
den2s=denin.*exp((Pac2s-Pac1in)./B);

aux=Cd;
ok=true;
for i=1:length(u1s)
    ok=true;
    
    for j=1:6

    dP1dt=real(-(N*Cd*sqrt(2*((Pac1s(i)/den1s(i))-(Pac2s(i)/den2s(i)))...
        /(1-(N^2)*(m^2)))*Sag-u1s(i)*So)*B/(H1*St+Ho*So-So*y1s(i))); 
    Qag=(B*So*u1s(i)-dP1dt*(H1*St+Ho*So-So*y1s(i)))/(B*N);
    Re=abs(Dag.*Qag/(visc*Sag));
    Cd=CD2mod2(Re,relLD,m);
    end
    Cds(i,1)=Cd;
    Qags(i,1)=Qag;

end


%%%===>     BAJADA DEL PISTON

%%%% Actualizacion de las condiciones iniciales

CI(1)=y1s(length(y1s));
CI(2)=u2s(length(u2s));
CI(3)=y2s(length(y2s));
CI(4)=Pac1s(length(Pac1s));
CI(5)=Pac2s(length(Pac2s));

%%%% Actualizacion de los parametros

%NOTA: Los parametros del 1 al 16 se mantienen
param(19)=PN2s(length(PN2s));
param(20)=CI(3);
param(21)=CI(4);
param(22)=CI(5);
param(23)=den1s(end);
param(24)=den2s(end);
param(25)=B;
param(26)=visc;
param(27)=Fc;
param(28)=Fs;
param(29)=vs;
param(30)=n;
param(31)=Fv;

[tb,yb]=ode23s(@(t2,x) sistodesbaj(t2,x,param),[0 tfin], CI);

y1b=yb(:,1);
u1b=A*W*sin(W*tb+pi);
u2b=yb(:,2);
y2b=yb(:,3);

%%%% Recalculo del resto de magnitudes para la bajada
VN2b=(St*(Ht-Eemb-y2b));
PN2b= PN2in.*((VN2in./VN2b).^(gamma));
Pac1b=yb(:,4);
Pac2b=yb(:,5);
den1b=den1s(end).*exp((Pac1b-CI(4))./B);
den2b=den2s(end).*exp((Pac2b-CI(5))./B);
clear dP2dt
for i=1:length(tb)
   if i>1 & i<length(tb)
        dP2dt(i,1)=(Pac2b(i+1)-Pac2b(i-1))/(tb(i+1)-tb(i-1));
   end
   if i==1
       dP2dt(1,1)=(Pac2b(2)-Pac2b(1))/(tb(2)-tb(1));
   end
   if i==length(tb)
       dP2dt(i,1)=(Pac2b(end)-Pac2b(end-1))/(tb(end)-tb(end-1));
   end
end

Qagb=-(B.*St*u2b+dP2dt.*St.*(y2b-Z))./(N*B); 


%%%% Se unen los vectores de subida y bajada
t=[ts' ts(end)+tb'];
y1=[y1s' y1b'];
u1=[u1s' u1b'];
y2=[y2s' y2b'];
u2=[u2s' u2b'];
PN2=[PN2s' PN2b'];
Pac1=[Pac1s' Pac1b'];
Pac2=[Pac2s' Pac2b'];
den1=[den1s' den1b'];
den2=[den2s' den2b'];
Qag=[Qags' -Qagb'];
Rec=Dag.*Qag/(visc*Sag);

%El siguiente fragmento de codigo se utiliza para que los vectores de las
%variables calculadas tengan la misma longitud
taux=t;
t=[0:0.0001:tfinal];
[taux,index]=unique(taux);
y2=interp1(taux,y2(index),t);
y1=interp1(taux,y1(index),t);
u2=interp1(taux,u2(index),t);
u1=interp1(taux,u1(index),t);
Pac1=interp1(taux,Pac1(index),t);
Pac2=interp1(taux,Pac2(index),t);
PN2=interp1(taux,PN2(index),t);
den1=interp1(taux,den1(index),t);
den2=interp1(taux,den2(index),t);
Qag=interp1(taux,Qag(index),t);
Rec=interp1(taux,Rec(index),t);


%Fuerza de fricción 
clear Ff1
    for j=1:length(u1)
        if u1(j)>=0
            Ff1(j) = Fc+(Fs-Fc).*exp(-(abs((u1(j)./vs))).^2)+Fv*u1(j);
            
        else
            Ff1(j) = -(Fc+(Fs-Fc).*exp(-(abs((u1(j)./vs))).^2))+Fv*u1(j);
            
        end
    end

    
%%%%%Calculo de la fuerza de amortiguación
%Una vez obtenida la presion P1 mediante el sistema, 
%es directo calcular esta fuerza

Famortcomp=(Pac1-Pac1(1))*So+Ff1;

