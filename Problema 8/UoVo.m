function dxdt  = UoVo(t,x,Uinf,Te1,Gamma_o,Variable)
% En esta función se introduce las derivadas de primer orden 
% de las variables que se desean calcular. 

% t:tiempo
% x: Variables que se integrarán
%    => x(1): Coordenada horizontal del vórtice
%    => x(2): Coordenada vertical del vortice
% dxdt: Derivadas de primer orden de cada una de las variables: 
%    => dxdt(1): Velocidad horizontal del vórtice
%    => dxdt(2): Velocidad vertical del vórtice
% Uinf: Velocidad no perturbada 
% Te1: Constante de tiempo
% Gamma_o: Circulación inicial si se considera variable.
%          Si la circulación es constante, se considerara este valor
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

