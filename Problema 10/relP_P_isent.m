function P2P1 = relP_P_isent( M1,M2,gamma )
%Relaci�n de presiones entre dos puntos si el flujo es isentr�pico
%M: N�mero de mach en el punto 
%gamma: �ndice de politrop�a
%P2P1: Relaci�n entre la presi�n aguas abajo (P2) y aguas arriba (P1)

P2P1=((1+(1/2)*M1^2*gamma-(1/2)*M1^2)^(gamma/(gamma-1)))/...
        ((1+(1/2)*M2^2*gamma-(1/2)*M2^2)^(gamma/(gamma-1)));

end

