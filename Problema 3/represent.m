function [ ] = represent( t,y1,u1,y2,u2,Pac1,Pac2,PN2,den1,den2, Qag, Re,Famort,Ff1,Do,difdiam)
% Representación gráfica
k=1;
TFuente=12; %Tamaño de la fuente
PointaxX=6; %Número de puntos  eje X
PointaxY=7; %Número de puntos  eje Y

if difdiam==false

p.fig=figure(k);
p.fig.Position=[1 93 757 635];

subplot(3,1,1)
p.dib(1)=plot(t,Pac1,'LineWidth',1.5);
hold on
grid on
p.dib(2)=plot(t,Pac2,'LineWidth',1.5);
p.dib(3)=plot(t,PN2,'LineWidth',1.5);

p.dib(1).Marker='^';
p.dib(2).Marker='o';
p.dib(3).Marker='*';

npuntos=25;
distpunt=fix(length(t)/25);
p.dib(1).MarkerIndices=(2:distpunt:length(t));
p.dib(2).MarkerIndices=(10:distpunt:length(t));
p.dib(3).MarkerIndices=(2:distpunt:length(t));

grid minor
axvect=axis;
ax=gca;
ax.XTick=(linspace(axvect(1),axvect(2),PointaxX));
ax.YTick=(linspace(axvect(3),axvect(4),5));
ax.FontSize=TFuente;
ax.YAxis.TickLabelFormat='%.3f';
ax.XAxis.TickLabelFormat='%.3f';

xlabel('Tiempo(s)')
ylabel('Presión (Pa)')
title('Presión cámara inferior (P1) e intermedia (P2)')
legend('P1','P2','PN2')

subplot(3,1,2)
plot(t,u1,'LineWidth',1.5)
hold on
grid on
grid minor

axvect=axis;
ax=gca;
ax.XTick=(linspace(axvect(1),axvect(2),PointaxX));
ax.YTick=(linspace(axvect(3),axvect(4),5));
ax.FontSize=TFuente;
ax.YAxis.TickLabelFormat='%.3f';
ax.XAxis.TickLabelFormat='%.3f';

xlabel('Tiempo (s)')
ylabel('Velocidad (m/s)')
title('Velocidad del pistón')

subplot(3,1,3)
plot(t,Qag,'LineWidth',1.5)
hold on; grid on; grid minor;

axvect=axis;
ax=gca;
ax.XTick=(linspace(axvect(1),axvect(2),PointaxX));
ax.YTick=(linspace(axvect(3),axvect(4),5));
ax.FontSize=TFuente;
ax.YAxis.TickLabelFormat='%.3f';
ax.XAxis.TickLabelFormat='%.3f';

xlabel('Tiempo (s)')
ylabel('Caudal (m^3/s)')
title('Caudal circulante a través de cada agujero')

k=k+1;
p=figure(k);
p.Position=[415 19 569 424];

plot(t,Famort,'LineWidth',1.8)
hold on; grid on; grid minor;

axvect=axis;
clear ax
ax=gca;
ax.XTick=(linspace(axvect(1),axvect(2),PointaxX));
ax.YTick=(linspace(axvect(3),axvect(4),PointaxY));
ax.FontSize=TFuente;
ax.YAxis.TickLabelFormat='%.1f';
ax.XAxis.TickLabelFormat='%.3f';

xlabel('Tiempo (s)')
ylabel('Fuerza de amortiguamiento (N)')

else
    p.fig=figure;
    p.fig.Position=[721 24 560 420];
    for i=1:length(Do)
        p.dib(i)=plot(t,Famort(i,:),'LineWidth',1.5);
        hold on
    end

    p.dib(1).Marker='*';
    p.dib(1).MarkerIndices=(1:25:length(t));
    p.dib(2).Marker='+';
    p.dib(2).MarkerIndices=(1:25:length(t));
    p.dib(3).Marker='x';
    p.dib(3).MarkerIndices=(1:25:length(t));
    p.dib(4).Marker='o';
    p.dib(4).MarkerIndices=(1:25:length(t));

    ax=gca;
    axvect=axis;
    ax.XTick=(linspace(axvect(1),axvect(2),6));
    ax.YTick=(linspace(axvect(3),axvect(4),11));
    ax.FontSize=12;
    ax.YAxis.TickLabelFormat='%.1f';
    ax.XAxis.TickLabelFormat='%.3f';

    grid on; grid minor;
    xlabel('Tiempo (s)')
    ylabel('Fuerza de amortiguamiento (N)')

    legend(['Do=' num2str(Do(1)) ' m'],['Do=' num2str(Do(2)) ' m'],...
           ['Do=' num2str(Do(3)) ' m'],['Do=' num2str(Do(4)) ' m'])

    p.fig=figure;
    p.fig.Position=[720 308 560 420];
    for i=1:length(Do)
        p.dib(i)=plot(t,Ff1(i,:),'LineWidth',1.5);
        hold on
    end
    vector(:,1)=(2:50:500);
    vector(:,2)=(36:50:500);
    vector(:,3)=(24:50:500);
    vector(:,4)=(12:50:500);

    p.dib(1).Marker='*';
    p.dib(1).MarkerIndices=vector(:,1);
    p.dib(2).Marker='+';
    p.dib(2).MarkerIndices=vector(:,2);
    p.dib(3).Marker='x';
    p.dib(3).MarkerIndices=vector(:,3);
    p.dib(4).Marker='o';
    p.dib(4).MarkerIndices=vector(:,4);

    ax=gca;
    axvect=axis;
    ax.XTick=(linspace(axvect(1),axvect(2),PointaxX));
    ax.YTick=(linspace(axvect(3),axvect(4),PointaxY));
    ax.FontSize=12;
    ax.YAxis.TickLabelFormat='%.1f';
    ax.XAxis.TickLabelFormat='%.3f';

    grid on; grid minor;
    xlabel('Tiempo (s)')
    ylabel('Fuerza de fricción (N)')

    legend(['Do=' num2str(Do(1)) ' m'],['Do=' num2str(Do(2)) ' m'],...
           ['Do=' num2str(Do(3)) ' m'],['Do=' num2str(Do(4)) ' m'])

end
end

