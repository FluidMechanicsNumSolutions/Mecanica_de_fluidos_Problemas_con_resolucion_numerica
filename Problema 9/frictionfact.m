function [f] = frictionfact(Re,rE)
%Permite obtener el factor de fricci�n  
% Re: N�mero de Reynolds
% rE: Rugosidad relativa
% f: Factor de fricci�n
options=optimoptions('fsolve','Display','off');
if Re<=2500 % Condici�n de fluido viscoso a bajas velocidades
    sol= 64/Re; %Seg�n Poiseulle 
    f=sol;
    
else if 4500<Re<10000 & rE==0  %Condici�n de flujo turbulento
                               %sin rugosidad (superficie suave)
    
    sol=0.316/(Re^(1/4));%Seg�n Blasius
    f=sol;
    end
    
end
if Re>10000 & rE==0 %Condici�n para flujo turbulento 
                    %completamente desarrollado sin rugosidad  
    x0=[0.1;0];
    sol=fsolve(@(x) f1(x,Re),x0,options);
    f=sol(1);
    return
end
if 4500<Re & Re<10000 & ne(rE,0) % Condici�n de flujo turbulento 
                                 % con rugosidad
    x0=[0.01;0];
    sol=fsolve(@(x) f2(x,Re,rE),x0,options);
    f=sol(1);
    return
end
if Re>10000 & ne(rE,0) % Condici�n para flujo turbulento 
                       % completamente desarrollado con rugosidad 
    x0=[0.1;0];
    sol=fsolve(@(x) f3(x,Re,rE),x0,options);
    f=sol(1);
    return
end

end

