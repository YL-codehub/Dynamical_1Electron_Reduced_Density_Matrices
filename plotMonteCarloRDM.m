  function RDMMap = plotMonteCarloRDM(solPs,atomicOrbitals)
%% This functions aims to plot DCPs performances : refined vs data ones, differences and anisotropies.
%
% --- INPUTS --- 
% - 'solPs' ((nOrbs x nOrbs) x N) : N concatenated population matrices from
% the same SDP refinement but with different noises on data.
% - 'atomicOrbitals' (structure) : non-orthogonalized atomic orbitals.
% ---------------
%
% --- OUTPUTS ---
% - 'RDMMap' : grid of the mean map.
% ---------------
%
% YL.

%% Compute maps
x=-5:0.1:5;
y=-3:0.1:3;
[~,path] = createMap(x,y);

atomicOrbitalsValuesOnPath = zeros(length(atomicOrbitals),size(path,1));

for i = 1:length(atomicOrbitals)
   atomicOrbitalsValuesOnPath(i,:) = atomicOrbitals{i}.eval(path); 
end

C = orthogonalisation(atomicOrbitals);

orthogAtomicOrbitalsValuesOnPath = C*atomicOrbitalsValuesOnPath;
n = length(atomicOrbitals);
N = fix(length(solPs)/n);
RDMMap_min = orthogAtomicOrbitalsValuesOnPath.'*solPs(:,1:n)*orthogAtomicOrbitalsValuesOnPath;
RDMMap_max = RDMMap_min;
RDMMap_mean = RDMMap_min./N;
RDMMap_sigma = zeros(n,n);

for i=(n+1):n:(length(solPs))
    % Uncomment/comment if you want maps with either minimum or maximum points :
    % RDMMap_mean = min(orthogAtomicOrbitalsValuesOnPath.'*solPs(:,i:(n-1+i))*orthogAtomicOrbitalsValuesOnPath,RDMMap_min);
    % RDMMap_mean = max(orthogAtomicOrbitalsValuesOnPath.'*solPs(:,i:(n-1+i))*orthogAtomicOrbitalsValuesOnPath,RDMMap_max);
    RDMMap_mean = RDMMap_mean + orthogAtomicOrbitalsValuesOnPath.'*solPs(:,i:(n-1+i))*orthogAtomicOrbitalsValuesOnPath./N;
end

[a,b] = size(RDMMap_mean);
MAPS = zeros(N,a,b);
for k= 1:N
    MAPS(k,:,:) = orthogAtomicOrbitalsValuesOnPath.'*solPs(:,1+(k-1)*n:k*n)*orthogAtomicOrbitalsValuesOnPath;
end

RDMMap_sigma = reshape(std(MAPS), [a,b]);
RDMMap_fluc = RDMMap_sigma;
RDMMap = RDMMap_fluc;


%% plot standard deviations at each point
valuesP = zeros(1,21);
for i=0:10
    valuesP(i+1)=0.001*2^i;
end
values = valuesP;
values = [-fliplr(valuesP),valuesP];
RDMP=RDMMap;
RDMP(RDMMap<=0) =0;
RDMM=RDMMap;
RDMM(RDMMap>=0) =0;
X=vecnorm(path-path(1,:),2,2);
figure;
contourf(X,X,RDMP,'LevelList',values,'LineWidth',0.1)
hold on
colormap(jet)
c = colorbar('southoutside');
c.Label.String = 'Density error ( per cube Angströms)';
set(gcf, 'Position',  [500, 500, 500, 500])
axis equal
xlabel('(a.u.)');
ylabel('(a.u.)');
set(gca,'xtick',0:ceil(X(end))-1)
set(gca,'ytick',0:ceil(X(end))-1)
set(gca,'ColorScale','log')
% colorbar();
plot(1,1,'o-','MarkerFaceColor','black','MarkerEdgeColor','black')
text(1+0.3,1+0.3, 'O','FontSize',18,'FontWeight','bold');
plot(3,3,'o-','MarkerFaceColor','black','MarkerEdgeColor','black')
text(3+0.3,3+0.3, 'C','FontSize',18,'FontWeight','bold');
plot(5,5,'o-','MarkerFaceColor','black','MarkerEdgeColor','black')
text(5+0.3,5+0.3, 'O','FontSize',18,'FontWeight','bold');

set(findall(gcf,'-property','FontSize'),'FontSize',15)

%% Plot mean RDM
hold off;
RDMMap = RDMMap_mean;
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

set(findall(gcf,'-property','FontSize'),'FontSize',20)
end

