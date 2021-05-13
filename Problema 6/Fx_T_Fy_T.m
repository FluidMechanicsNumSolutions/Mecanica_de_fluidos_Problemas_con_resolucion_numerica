function [Fx_T,Fy_T] = Fx_T_Fy_T(delta,e,mu,w,R,h)
%Función que calcula la fuerza horizontal y vertical sobre el eje
% delta: Coordenada angular del eje respecto el centro
% e: excentricidad
% mu: viscosidad dinámica
% w: velocdidad angular
% R: radio del eje
% h: Diferencia entre el radio del estator y del eje
% Fx_T: Fuerza horizontal sobre el eje
% Fy_T: Fuerza vertical sobre el eje

psi=h/R;
epsilon=e/h;

% Integrales I1, I2, I3, I4, I5, I6, I7

I(1)=2*pi/sqrt(1-epsilon^2);
I(2)=2*pi/((1-epsilon^2)^(3/2));
I(3)=pi*(2+epsilon^2)/((1-epsilon^2)^(5/2));
I(4)=(I(2)-I(1))/epsilon; 
I(5)=(I(3)-I(2))/epsilon;
I(6)=0; 
I(7)=0;

% Fuerza horizontal y vertical sobre el eje

Fx_cil=(6*mu*w*R/psi^2)*(I(6)*I(3)-I(7)*I(2))/I(3);
Fy_cil=(6*mu*w*R/psi^2)*(I(2)*I(5)-I(3)*I(4))/I(3);

if isnan(Fy_cil)
    Fy_cil=0;
end

Fx_T=-Fy_cil*sin(delta)-Fx_cil*cos(delta);
Fy_T=Fy_cil*cos(delta)-Fx_cil*sin(delta);


end

