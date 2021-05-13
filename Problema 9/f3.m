function fc = f3(f,Re,rE)
% Ecuación definida por Karman-Prandtl 
    x=f(1);
    fc(1)=2*log10(1/(2*rE))+1.74-1/sqrt(x);
    fc(2)=0;
end
