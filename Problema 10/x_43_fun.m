function [M,P,T,F,Fx,Fy,drag,lift,power] = x_43_fun(X,Mi,alpha,gamma,Po,To,R)

% Esta función calcula las variables de salida que se piden en el problema
% apartados 1,2,3. El objeto de su definición es poder utilizar el mismo
% algoritmo para calcular los apartados 1,2,3 y los diferentes casos del
% apartado 4. 

%   X: Coordenadas horizontales X1,X2,X3,X4,X5
%   Mi: Número de Mach en la entrada 
%   alpha: Ángulos de inclinación
%   gamma: Constante politrópica
%   Po: Presión a una altura determinada
%   To: Temperatura a una altura determinada

%   M: Numeros de Mach en cada una de las superficies
%   P: Presión estática en cada superfície
%   T: Temperatura en cada superficie
%   F: Fuerzas generadas sobre cada superficie
%   Fx: Fuerzas de arrastre asociadas a cada una de las superficies
%   Fy: Fuerzas de sustentación asociadas a cada una de las superficies
%   drag: Fuerza de arrastre total
%   lift: Fuerza de sustentación total

    %Inicialización de variables
    M=zeros(5,1);
    P=zeros(5,1);
    T=zeros(5,1);
    F=zeros(5,1);
    Fx=zeros(5,1);
    Fy=zeros(5,1);

    % Determinación de las características del fluido en la parte superior:
    % => Superficie 1
    epsilon(1)=f_epsilon(alpha(1),gamma,Mi);
    Mn_sup=Mi*sin(epsilon(1));

    Mn(1)=MNumNormShock(Mn_sup,gamma);
    M(1)=Mn(1)/sin(epsilon(1)-alpha(1));
    P(1)=Po*r2PNormShock(Mn_sup,gamma); 
    T(1)=To*r2TNormShock(Mn_sup,gamma);
    
    % Determinación de las características del fluido en la parte inferior:
    % => Superficie 2

    epsilon(2)=f_epsilon(alpha(2),gamma,Mi);
    Mn_inf=Mi*sin(epsilon(2));

    Mn(2)=MNumNormShock(Mn_inf,gamma);
    M(2)=Mn(2)/sin(epsilon(2)-alpha(2));

    P(2)=Po*r2PNormShock(Mn_inf,gamma); 
    T(2)=To*r2TNormShock(Mn_inf,gamma);

    % => Superficie 3

    epsilon(3)=f_epsilon(alpha(3)-alpha(2),gamma,M(2));
    Mn_23=M(2)*sin(epsilon(3));

    Mn(3)=MNumNormShock(Mn_23,gamma);
    M(3)=Mn(3)/sin(epsilon(3)-(alpha(3)-alpha(2)));

    P(3)=P(2)*r2PNormShock(Mn_23,gamma); 
    T(3)=T(2)*r2TNormShock(Mn_23,gamma);

    % => Superficie 4

    eta34_inf=f_eta(gamma,M(3),"eta"); 
    eta34_inf_2=eta34_inf+alpha(3); 

    M(4)=f_eta(gamma,M(3),"M",eta34_inf_2);
    
    P(4)=P(3)*relP_P_isent(M(3),M(4),gamma);
    T(4)=T(3)*relT_T_isent(M(3),M(4),gamma);

    % => Superficie 5

    eta45_inf=f_eta(gamma,M(4),"eta"); 
    eta45_inf_1=eta45_inf+alpha(3); 

    M(5)=f_eta(gamma,M(4),"M",eta45_inf_1);

    P(5)=P(4)*relP_P_isent(M(4),M(5),gamma);
    T(5)=T(4)*relT_T_isent(M(4),M(5),gamma);

    % Fuerzas sobre la aeronave 

    for i=1:5
       if i<4
        F(i)=P(i)*X(i)/cos(alpha(i)); 
        Fx(i)=F(i)*sin(alpha(i));
       if i==1
           Fy(i)=-F(i)*cos(alpha(i));
       else
           Fy(i)=F(i)*cos(alpha(i));
       end
       else
           if i==4
               F(4)=P(4)*X(4);
               Fx(4)=0;
               Fy(4)=F(4);
           else
               F(5)=P(5)*X(5)/cos(alpha(4));
               Fx(5)=-F(5)*sin(alpha(4));
               Fy(5)=F(5)*cos(alpha(4));
           end
       end

    end

    drag=sum(Fx);
    lift=sum(Fy); 
    
    % Potencia necesaria de los motores
    V=Mi*sqrt(gamma*R*To);
    power=drag*V;
    
end




