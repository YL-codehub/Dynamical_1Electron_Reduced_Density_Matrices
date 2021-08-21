  function RDMMap = plotRDM_thermaldef(P1,atomicOrbitals1,P2,atomicOrbitals2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x=-5:0.1:5;
y=-3:0.1:3;
[~,path] = createMap(x,y);

atomicOrbitalsValuesOnPath = zeros(length(atomicOrbitals1),size(path,1));
atomicOrbitalsValuesOnPath2 = zeros(length(atomicOrbitals2),size(path,1));

for i = 1:length(atomicOrbitals1)
   atomicOrbitalsValuesOnPath(i,:) = atomicOrbitals1{i}.eval(path); 
   atomicOrbitalsValuesOnPath2(i,:) = atomicOrbitals2{i}.eval(path); 
end

C = orthogonalisation(atomicOrbitals1);
C2 = orthogonalisation(atomicOrbitals2);

orthogAtomicOrbitalsValuesOnPath = C*atomicOrbitalsValuesOnPath;
orthogAtomicOrbitalsValuesOnPath2 = C2*atomicOrbitalsValuesOnPath2;

RDMMap = orthogAtomicOrbitalsValuesOnPath.'*P1*orthogAtomicOrbitalsValuesOnPath-orthogAtomicOrbitalsValuesOnPath2.'*P2*orthogAtomicOrbitalsValuesOnPath2;



%% plot
valuesP = zeros(1,21);
for i=0:20
    valuesP(i+1)=0.001*2^i;
end
values = [-fliplr(valuesP),valuesP];
RDMP=RDMMap;
RDMP(RDMMap<=0) =0;
RDMM=RDMMap;
RDMM(RDMMap>=0) =0;
X=vecnorm(path-path(1,:),2,2);
figure;
contour(X,X,RDMP,'LevelList',values,'LineColor','b','LineWidth',0.9);
hold on;
contour(X,X,RDMM,'LevelList',values,'LineColor','r','LineStyle','--','LineWidth',0.9);
hold on;
contour(X,X,RDMP,'LevelList',[0],'LineColor','g','LineStyle','--','LineWidth',0.9);
hold on;
set(gcf, 'Position',  [500, 500, 500, 500]);
axis equal;
xlabel('(a.u.)');
ylabel('(a.u.)');
set(gca,'xtick',0:ceil(X(end))-1);
set(gca,'ytick',0:ceil(X(end))-1);

plot(1,1,'o-','MarkerFaceColor','black','MarkerEdgeColor','black');
text(1+0.3,1+0.3, 'O','FontSize',18,'FontWeight','bold');
plot(3,3,'o-','MarkerFaceColor','black','MarkerEdgeColor','black');
text(3+0.3,3+0.3, 'C','FontSize',18,'FontWeight','bold');
plot(5,5,'o-','MarkerFaceColor','black','MarkerEdgeColor','black');
text(5+0.3,5+0.3, 'O','FontSize',18,'FontWeight','bold');
set(findall(gcf,'-property','FontSize'),'FontSize',20)
end

