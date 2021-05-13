function [ Msol ] = IsentMachNum(relA1A2,M1,gamma)
% Calcula el numero de Mach en un punto sabiendo la relación 
%de diámetros relA1A2=A2/A1, el número de Mach (M1) en uno de  
%los puntos y el índice de politropía gamma

Mach=[0.01:0.00001:10];
k=1;
Msol=0;
for i=1:length(Mach)
    M2=Mach(i);
    fun=[M1*((1+((gamma-1)*(1/2))*M2^2)/(1+((gamma-1)*(1/2))*...
        M1^2))^((gamma+1)/(2*(gamma-1)))/M2]-relA1A2;
    ec(i)=fun;
    if i==1
        funant=fun;
    else
        if ec(i)*ec(i-1)<0
            Msol(k,1)=M2;
            k=k+1;
        end
    end
end

end

