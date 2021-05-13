function P2P1 = relP_P_isent( M1,M2,gamma )
%Relación de presiones entre dos puntos si el flujo es isentrópico
%M: Número de mach en el punto 
%gamma: Índice de politropía
%P2P1: Relación entre la presión aguas abajo (P2) y aguas arriba (P1)

P2P1=((1+(1/2)*M1^2*gamma-(1/2)*M1^2)^(gamma/(gamma-1)))/...
        ((1+(1/2)*M2^2*gamma-(1/2)*M2^2)^(gamma/(gamma-1)));

end

