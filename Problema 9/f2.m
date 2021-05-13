function fc = f2(f,Re,rE)
% Ecuación definida por Colebrook
x=f(1);
fc(1)=1.14-2*log10(rE+9.35/(Re*sqrt(x)))-1/sqrt(x);
fc(2)=0;
end


