
function [du] = f2(t2,p)
% En esta función se introduce la derivada de primer orden de la presión
% sin considerar las fugas. 
% t2: tiempo
% p: presión
% du: derivada de la presión respecto del tiempo (sin fugas) 
 
global Rp Di Sp alpha Ao Vo l1 l2 l4 l6 l8 l10 l3 l5 l7 l9...
       l11 hs1 hs3 hs2 h1 h10 h2 h11 r1 r2 r3 r4 Pd w Cd B mu rho Ptank

du=B*(-Cd*Ao*sqrt(2)*sqrt((p-Pd)/rho)+Sp*Rp*tan(alpha)*w*sin(w*t2))/...
    (Vo+Sp*Rp*tan(alpha)*(cos(w*t2)+1));
end


