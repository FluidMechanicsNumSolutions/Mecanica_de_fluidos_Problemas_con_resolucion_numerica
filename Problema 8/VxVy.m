function [Vx,Vy,V,r] = VxVy(U,Gamma,X,Y,Xo,Yo)
% En esta función se calcula las velocidades del fluido (Vx, Vy) 
% en las direcciones X,Y teniendo en cuenta la posición del vortice Xo,Yo

% U: Velocidad no perturbada
% (Xo, Yo): Posición del vórtice
% Gamma: Circulación asociada al vórtice
% X,Y: Coordenadas de la posición a calcular
% r: radio del vórtice

Vx=U+Gamma/(2*pi)*((Y+Yo)./((X-Xo).^2+(Y+Yo).^2)-(Y-Yo)./...
    ((X-Xo).^2+(Y-Yo).^2)+(Y-Yo)./((X+Xo).^2+(Y-Yo).^2)-...
    (Y+Yo)./((X + Xo).^2+(Y+Yo).^2));

Vy=-Gamma/(2*pi)*((X-Xo)./((X-Xo).^2+(Y+Yo).^2)-(X-Xo)./...
    ((X-Xo).^2+(Y-Yo).^2)+(X+Xo)./((X+Xo).^2+(Y-Yo).^2)-...
    (X+Xo)./((X+Xo).^2+(Y+Yo).^2));

V=sqrt(Vx.^2+Vy.^2);

ind=find(V>25,1);
if length(V(:,1))>1 && length(V(1,:))>1
    r=sqrt((X(ind)-Xo)^2+(Y(ind)-Yo)^2);
    for i=1:length(V(:,1))
        for j=1:length(V(1,:))
            if ((X(i,j)-Xo)^2+(Y(i,j)-Yo)^2)<=(r)^2
               Vx(i,j)=NaN; 
               Vy(i,j)=NaN; 
               V(i,j)=NaN; 
            end
        end
    end
end
end

