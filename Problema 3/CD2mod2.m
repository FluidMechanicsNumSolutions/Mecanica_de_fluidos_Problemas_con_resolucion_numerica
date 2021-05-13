function Cd = CD2mod2(Re,relLD,m)
% Calculo del coeficiente de descarga. Fuente: Articulo 
% 'Discharge Coefficientos for incompressible non-cavitating flow through
% long orifices' Ecuaciones (4), (7) y (12) del articulo. 

%Utilizando metodo de regula-falsi 
diff=0.0001;
nodiff=true;
b=1;
a=0.5;
 fa=-1./a+1./(.827-0.85e-2*relLD)+20.*a.*(1+2.25.*relLD)...
     ./(sqrt(-m^2+1).*Re)-0.5e-2.*relLD...
     ./((1+7.5.*log(0.15e-3.*sqrt(-m^2+1).*Re./a).^2));
 fb=-1./b+1./(.827-0.85e-2*relLD)+20.*b.*(1+2.25.*relLD)...
     ./(sqrt(-m^2+1).*Re)-0.5e-2.*relLD...
     ./((1+7.5.*log(0.15e-3.*sqrt(-m^2+1).*Re./b).^2));
 p=(a*fb-b*fa)/(fb-fa);
 fp=-1./p+1./(.827-0.85e-2*relLD)+20.*p.*(1+2.25.*relLD)...
     ./(sqrt(-m^2+1).*Re)-0.5e-2.*relLD...
     ./((1+7.5.*log(0.15e-3.*sqrt(-m^2+1).*Re./p).^2));
if Re<1
    Re=1;
    fa=-1./a+1./(.827-0.85e-2*relLD)+20.*a.*(1+2.25.*relLD)...
        ./(sqrt(-m^2+1).*Re)-0.5e-2.*relLD...
        ./((1+7.5.*log(0.15e-3.*sqrt(-m^2+1).*Re./a).^2));
    fb=-1./b+1./(.827-0.85e-2*relLD)+20.*b.*(1+2.25.*relLD)...
        ./(sqrt(-m^2+1).*Re)-0.5e-2.*relLD...
        ./((1+7.5.*log(0.15e-3.*sqrt(-m^2+1).*Re./b).^2));
    p=(a*fb-b*fa)/(fb-fa);
    fp=0;
    Cd=p;
else
if Re<1e4 & Re>1

if abs(fp)<diff
    nodiff=false;
    Cd=p;
end  
while nodiff  
       
       if fa*fp<0
           b=p;
       else
           a=p;
       end
       fa=-1./a+1./(.827-0.85e-2*relLD)+20.*a.*(1+2.25.*relLD)...
           ./(sqrt(-m^2+1).*Re)-0.5e-2.*relLD...
           ./((1+7.5.*log(0.15e-3.*sqrt(-m^2+1).*Re./a).^2));
       fb=-1./b+1./(.827-0.85e-2*relLD)+20.*b.*(1+2.25.*relLD)...
           ./(sqrt(-m^2+1).*Re)-0.5e-2.*relLD...
           ./((1+7.5.*log(0.15e-3.*sqrt(-m^2+1).*Re./b).^2));
       p=(a*fb-b*fa)/(fb-fa);
       fp=-1./p+1./(.827-0.85e-2*relLD)+20.*p.*(1+2.25.*relLD)...
           ./(sqrt(-m^2+1).*Re)-0.5e-2.*relLD...
           ./((1+7.5.*log(0.15e-3.*sqrt(-m^2+1).*Re./p).^2));
       if abs(fp)<diff
           nodiff=false;
           Cd=p;
       end
    
end
else
Cd=0.827-0.85e-2*relLD;
end
end
end

