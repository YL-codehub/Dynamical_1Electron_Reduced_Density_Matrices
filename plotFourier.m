function [outputArg1,outputArg2] = plotFourier(map,SF,volume)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


x = map.points(:,1);
y = map.points(:,2);
z = map.points(:,3);
Qx = (SF.Q(:,1)).';
Qy = (SF.Q(:,2)).';
Qz = (SF.Q(:,3)).';

DENSList = exp(1i*(x.*Qx + y.*Qy + z.*Qz))*(SF.value(:,1)+1i*SF.value(:,2))/volume; 

DENSMap = reshape(DENSList,map.N,[]).';

%% plot
for i=0:20
    values(i+1)=0.01*2^(i);
end
values = [-fliplr(values),values];
DENSD = real(DENSMap);
DENSDP=DENSD;
DENSDP(DENSD<=0) =0;
DENSDM=DENSD;
DENSDM(DENSD>=0) =0;

Y=vecnorm(map.points(1:map.N,:)-map.points(1,:),2,2);
X=vecnorm(map.points(1:map.N:end,:)-map.points(1,:),2,2);
figure;
contour(Y(1:end-1),X,DENSDP(:,1:end-1),'LevelList',values,'LineColor','b','LineWidth',1.1);
hold on;
contour(Y(1:end-1),X,DENSDM(:,1:end-1),'LevelList',values,'LineColor','r','LineStyle','--','LineWidth',1.1);
hold on;
axis equal;
xlabel('(a.u.)');
ylabel('(a.u.)');
set(gca,'xtick',0:ceil(Y(end))-1);
set(gca,'ytick',0:ceil(X(end))-1);
set(gcf, 'Position',  [500, 500, 500, 500]);

plot(3,2,'o-','MarkerFaceColor','black','MarkerEdgeColor','black');
text(3+0.3,2+0.3, 'O','FontSize',18,'FontWeight','bold');
plot(3,2+2,'o-','MarkerFaceColor','black','MarkerEdgeColor','black');
text(3+0.3,2+2+0.3, 'C','FontSize',18,'FontWeight','bold');
plot(3,2+4,'o-','MarkerFaceColor','black','MarkerEdgeColor','black');
text(3+0.3,2+4+0.3, 'O','FontSize',18,'FontWeight','bold');

set(findall(gcf,'-property','FontSize'),'FontSize',20);
end

