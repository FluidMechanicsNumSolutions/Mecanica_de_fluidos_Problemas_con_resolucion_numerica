function T2T1 = relT_T_isent( M1,M2,gamma )
%Relaci�n de temperaturas entre dos puntos si el flujo es isentr�pico
%M: N�mero de mach en el punto 
%gamma: �ndice de politrop�a
%T2T1: Relaci�n entre la temperatura aguas abajo (T2) y aguas arriba (T1)

T2T1=((1+(1/2)*M1^2*gamma-(1/2)*M1^2))/...
        ((1+(1/2)*M2^2*gamma-(1/2)*M2^2));

end

