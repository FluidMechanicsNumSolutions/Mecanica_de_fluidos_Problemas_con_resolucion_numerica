function dxdt  = UoVo(t,x,Uinf,Te1,Gamma_o,Variable)
% En esta funci�n se introduce las derivadas de primer orden 
% de las variables que se desean calcular. 

% t:tiempo
% x: Variables que se integrar�n
%    => x(1): Coordenada horizontal del v�rtice
%    => x(2): Coordenada vertical del vortice
% dxdt: Derivadas de primer orden de cada una de las variables: 
%    => dxdt(1): Velocidad horizontal del v�rtice
%    => dxdt(2): Velocidad vertical del v�rtice
% Uinf: Velocidad no perturbada 
% Te1: Constante de tiempo
% Gamma_o: Circulaci�n inicial si se considera variable.
%          Si la circulaci�n es constante, se considerara este valor
%          siempre
% Variable: Variable de tipo booleana. 
%       ==> Si la circulacion es variable: Variable = true
%       ==> Si la circulacion es constante: Variable = false

    x_vort=x(1);
    y_vort=x(2);
    if Variable
        Gamma=Gamma_o*exp(-t/(Te1));
    else
        Gamma=Gamma_o;
    end
    U=Uinf+(Gamma/(4*pi))*(x_vort^2)/(y_vort*(x_vort^2 + y_vort^2));
    V=-(Gamma/(4*pi))*y_vort^2/(x_vort*(x_vort^2 + y_vort^2));
    dxdt=[U;V];

end

