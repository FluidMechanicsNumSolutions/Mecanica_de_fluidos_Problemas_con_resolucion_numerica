function [du] = f3_down(t,p)
% En esta funci?n se introduce la derivada de primer orden de la presi?n
% considerando las fugas. (Durante el descenso)
% t: tiempo
% p: presi?n
% du: derivada de la presi?n respecto del tiempo (con fugas) 

global Rp Di Sp alpha Ao Vo l1 l2 l4 l6 l8 l10 l3 l5 l7 l9...
       l11 hs1 hs3 hs2 h1 h10 h2 h11 r1 r2 r3 r4...
       Pd w Cd B mu rho Ptank
   
QL1=pi*Di*((h1*Rp*tan(alpha)*sin(w*t)*w)*(1/2))-(1/12)*pi*Di*...
    (Ptank-p-(6*Rp*tan(alpha)*(sin(w*t)))*w*mu*(h10-h1)*...
    (l2+l4+l6+l8+l10)/h10^3)/(mu*(l1+l2+l3+l4+l5+l6+l7+l8+l9+...
    l10+l11-(0.195*10^(-1))*(1/2)-Rp*tan(alpha)*cos(w*t))/...
    (h11^3)+mu*(1/h2^3-1/h1^3)*(l2+l4+l6+l8+l10));
QL2=pi*(p-Ptank)/((6*mu)*(log(r2/r1)/hs1^3+log(r3/r2)/...
    hs2^3+log(r4/r3)/hs3^3));
du=B*(+Cd*Ao*sqrt(2)*sqrt((Pd-p)/rho)-Sp*Rp*tan(alpha)...
    *w*sin(w*t)-QL1-QL2)/(Vo+Sp*Rp*tan(alpha)*(1-cos(w*t)));

end

