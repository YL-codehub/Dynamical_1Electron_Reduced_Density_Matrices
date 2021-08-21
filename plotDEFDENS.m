function DENSMap = plotDEFDENS(map,P,atomicOrbitals,DENSPato)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
atomicOrbitalsValuesOnMap = zeros(length(atomicOrbitals),size(map.points,1));

for i = 1:length(atomicOrbitals)
   atomicOrbitalsValuesOnMap(i,:) = atomicOrbitals{i}.eval(map.points); 
end

C = orthogonalisation(atomicOrbitals);

orthogAtomicOrbitalsValuesOnMap = C*atomicOrbitalsValuesOnMap;
DENSList = diag(orthogAtomicOrbitalsValuesOnMap.'*P*orthogAtomicOrbitalsValuesOnMap);
DENSMap = reshape(DENSList,map.N,[]).';


%% plot
values = zeros(1,20);
for i=1:20
    values(i)=0.01*2^(i);
end
values = [-fliplr(values),-0.01,0.01,values];

DENSD = DENSMap-DENSPato;
DENSDP=DENSD;
DENSDP(DENSD<=0) =0;
DENSDM=DENSD;
DENSDM(DENSD>=0) =0;

figure;
Y=vecnorm(map.points(1:map.N,:)-map.points(1,:),2,2);
X=vecnorm(map.points(1:map.N:end,:)-map.points(1,:),2,2);
contour(Y(1:end-1),X,DENSDP(:,1:end-1),'LevelList',values,'LineColor','b','LineWidth',1.1)
hold on
contour(Y(1:end-1),X,DENSDM(:,1:end-1),'LevelList',values,'LineColor','r','LineStyle','--','LineWidth',1.1)
hold on
axis equal
xlabel('(a.u.)');
ylabel('(a.u.)');
set(gca,'xtick',0:ceil(Y(end))-1)
set(gca,'ytick',0:ceil(X(end))-1)
set(gcf, 'Position',  [500, 500, 500, 500])

plot(3,2,'o-','MarkerFaceColor','black','MarkerEdgeColor','black')
text(3+0.3,2+0.3, 'O','FontSize',18,'FontWeight','bold');
plot(3,2+2,'o-','MarkerFaceColor','black','MarkerEdgeColor','black')
text(3+0.3,2+2+0.3, 'C','FontSize',18,'FontWeight','bold');
plot(3,2+4,'o-','MarkerFaceColor','black','MarkerEdgeColor','black')
text(3+0.3,2+4+0.3, 'O','FontSize',18,'FontWeight','bold');

set(findall(gcf,'-property','FontSize'),'FontSize',20)

end
