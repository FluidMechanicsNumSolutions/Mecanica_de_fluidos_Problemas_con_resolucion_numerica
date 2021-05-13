function [ I ] = dI( alpha,params)
% Funci�n a integrar respecto del �ngulo
% alpha: �ngulo
% params: par�metros necesarios definidos en el script principal
% I: Expresi�n a integrar

Rr=params(1); R=params(2); H=params(3);
ga=params(4);

I=1./((2*R+H)*sin(ga+alpha)+(-2*H-2*R)*sin(ga)-2*Rr);

end

