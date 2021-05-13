function [r14,M1,M3,M4,P1,P2,P3,P4,T1,T2,T3,T4] = Subsonico(M2)
% En caso SUBSÓNICO, dado el número de Mach en 2,   
% se calcula el resto de magnitudes
% M1,M2,M3,M4: Número de Mach en 1,2,3,4 respectivamente
% P1,P2,P3,P4: Presión estática en 1,2,3,4 respectivamente [Pa]
% T1,T2,T3,T4: Temperatura estática 1,2,3,4 respectivamente [K]
% r14: relación de presión entre 1 y 4 

global L12 L34 D1 D3 Sup1 Sup3 rE12 rE34 R Po1 To1 gamma viscosdin

%%%% Cálculo del número de Mach en 1
% De la misma manera que se realizó al principio del programa en el script
% principal, se realiza un bucle para determinar M1 junto con los valores
% correctos del número de Reynolds y el factor de fricción asociado

M1=0.5; %Suposición inicial
Re1=massflow(M1,Po1,To1,Sup1,R,gamma)*D1/(Sup1*viscosdin);
i=0;
f1=frictionfact(Re1,rE12);
M1new=FannoMachNum(f1*L12/D1+FannoIndvar(M2,gamma),gamma);
if length(M1new)>1
    M1new=M1new(1);
end
while abs(M1new-M1)>0.0001 
    M1=M1new;
    Re1=massflow(M1,Po1,To1,Sup1,R,gamma)*D1/(Sup1*viscosdin);
    f1=frictionfact(Re1,rE12);
    M1new=FannoMachNum(f1*L12/D1+FannoIndvar(M2,gamma),gamma);
    i=i+1;
end
if length(M1new)>1
    M1new=M1new(1);
end
M1=M1new;
T1=To1/relToT(M1,gamma);
P1=Po1/relPoP(M1,gamma);

%%%% Relacion entre las magnitudes en el punto 1 y 2 (Flujo de Fanno)
P2=(P1/Fannorel2P(M1,gamma))*(Fannorel2P(M2,gamma));
Po2=(Po1/Fannorel2Po(M1,gamma))*(Fannorel2Po(M2,gamma));
T2=(T1/Fannorel2T(M1,gamma))*Fannorel2T(M2,gamma);
To2=To1;

%%%% Cálculo del número de Mach en 3
relA2A3=Sup3/Sup1; 
M3c=IsentMachNum(relA2A3,M2,gamma);
% Se obtienen dos soluciones para M3. Se descartara aquella cuyo número de
% Mach sea superior a la unidad dado que en este caso se estudia el caso
% subsónico

Po3=Po2;
To3=To2;
To4=To3;

for i=1:length(M3c) %M3 debe ser forzosamente subsónico
    if M3c(i)<1
        M3=M3c(i);
    end
end
T3=To3/relToT(M3,gamma);
P3=Po3/relPoP(M3,gamma);

%%%% Cálculo del numero de Mach en 4 

% M4=[];
fLD3=FannoIndvar(M3,gamma);
Re3=massflow(M3,Po3,To3,Sup3,R,gamma)*D3/(Sup3*viscosdin);
f34=frictionfact(Re3,rE34);
varindep=(fLD3-(f34*L34/D3));
M4c=FannoMachNum(varindep,gamma);
% ¡Atención! En el caso que no sea posible obtener un número de Mach en
% 4 para uno de los números de Mach calculados en 3, el resultado de la
% anterior función sera NaN. 

%%%% Evaluación de relaciones de presion entre 1 y 4 (Po1/P4)

M4=NaN;
for i=1:length(M4)
    if M4c(i)<1 & ne(isnan(M4c(i)),true)
        M4=M4c(i);
    end
end
r14=(Po1/Po2)*(Po2/Po3)*Fannorel2Po(M3,gamma)*...
    (1/Fannorel2Po(M4,gamma))*relPoP(M4,gamma);
T4=To4/relToT(M4,gamma);
P4=Po1/r14;


end

