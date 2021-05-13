function [ Msol ] = FannoMachNum( varind,gamma )
% Se obtiene el número de Mach de un punto tomando como referencia 
% otro punto en condición critica. 
% Esta función es válida para flujo de Fanno 

% varind= 4*Cf*Ax/D
% Msol: Número de Mach resultante
ec=[];
Mach=(0.01:0.00001:10);
Msol=0;
k=1;
if varind>0
for i=1:length(Mach)
    M=Mach(i);
    fun=((-M^2+1)/(gamma*M^2)+(gamma+1)*log((gamma+1)*M^2/...
        (2*(1+((gamma-1)*(1/2))*M^2)))/(2*gamma))-varind;
    ec(i)=fun;
    if i==1
        funant=fun;
    else
        if ec(i)*ec(i-1)<0
            Msol(k,1)=M;
            k=k+1;
        end
    end
end
if Msol==0
    Msol=NaN;
end
else
    if varind==0
        Msol=1;
    else
        Msol=NaN;
    end
end
end

