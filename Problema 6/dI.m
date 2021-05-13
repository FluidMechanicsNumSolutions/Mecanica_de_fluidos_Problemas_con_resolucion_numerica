function [ I ] = dI( alpha,params)
% Función a integrar respecto del ángulo
% alpha: ángulo
% params: parámetros necesarios definidos en el script principal
% I: Expresión a integrar

Rr=params(1); R=params(2); H=params(3);
ga=params(4);

I=1./((2*R+H)*sin(ga+alpha)+(-2*H-2*R)*sin(ga)-2*Rr);

end

