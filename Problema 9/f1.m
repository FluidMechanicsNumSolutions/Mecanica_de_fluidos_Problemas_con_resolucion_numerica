function [fc] = f1(f,Re)
% Ecuación definida por Karman-Prandtl 
fc(1)=2.*log10(Re.*(f(1).^(1/2)))-0.8-1./(f(1).^(1/2));
fc(2)=0;
end

