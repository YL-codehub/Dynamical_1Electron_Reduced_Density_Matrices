function RDMMap = plotRDM(path,P,atomicOrbitals)
%   function RDMMap = plotRDM(path,P,atomicOrbitals)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
atomicOrbitalsValuesOnPath = zeros(length(atomicOrbitals),size(path,1));

for i = 1:length(atomicOrbitals)
   atomicOrbitalsValuesOnPath(i,:) = atomicOrbitals{i}.eval(path); 
end

C = orthogonalisation(atomicOrbitals);

orthogAtomicOrbitalsValuesOnPath = C*atomicOrbitalsValuesOnPath;

RDMMap = orthogAtomicOrbitalsValuesOnPath.'*P*orthogAtomicOrbitalsValuesOnPath;

%% plot
valuesP = zeros(1,21);
for i=0:20
    valuesP(i+1)=0.01*2^i;
end
values = [-fliplr(valuesP),valuesP];
RDMP=RDMMap;
RDMP(RDMMap<=0) =0;
RDMM=RDMMap;
RDMM(RDMMap>=0) =0;
X=vecnorm(path-path(1,:),2,2);
figure;
contour(X,X,RDMP,'LevelList',values,'LineColor','b','LineWidth',0.9)
hold on
contour(X,X,RDMM,'LevelList',values,'LineColor','r','LineStyle','--','LineWidth',0.9)
hold on
contour(X,X,RDMP,'LevelList',[0],'LineColor','g','LineStyle','--','LineWidth',0.9)
hold on
set(gcf, 'Position',  [500, 500, 500, 500])
axis equal
xlabel('(a.u.)');
ylabel('(a.u.)');
set(gca,'xtick',0:ceil(X(end))-1)
set(gca,'ytick',0:ceil(X(end))-1)

plot(1,1,'o-','MarkerFaceColor','black','MarkerEdgeColor','black')
text(1+0.3,1+0.3, 'O','FontSize',18,'FontWeight','bold');
plot(3,3,'o-','MarkerFaceColor','black','MarkerEdgeColor','black')
text(3+0.3,3+0.3, 'C','FontSize',18,'FontWeight','bold');
plot(5,5,'o-','MarkerFaceColor','black','MarkerEdgeColor','black')
text(5+0.3,5+0.3, 'O','FontSize',18,'FontWeight','bold');
%title(['$$\chi^2/ndf = $$'],'interpreter','latex');

end

