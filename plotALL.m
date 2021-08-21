%function [SFpredicted,CPpredicted,RDMMap] = plotALL(P,atomicOrbitals,R,T,SF,CP,volume,pathPato, B)
function [rdm,dens] = plotALL(P,atomicOrbitals,R,T,SF,CP,volume,pathPato, B)
%% This functions aims to plot DCPs performances : refined vs data ones, differences and anisotropies.
%
% --- INPUTS --- 
% - 'P' (nOrbs x nOrbs) : Population matrix from refinements (still
% orthogonalized)
% - 'atomicOrbitals' (structure) : non-orthogonalized atomic orbitals.
% - 'R' (nMolecules x 3x3) : rotation matrices to go from reference molecules
% to all in the cell.
% - 'T' (nMolecules x 3) : DCPs data points.
% - 'SF' (structure) : SFs data points.
% - 'CP' (structure) : DCPs data points.
% - 'volume' (float) : volume of the primitive cell given by readCoords
% function.
% - 'pathPato' (float) : .dat pseudo-atomic density map file for
% deformation map.
% - 'B' (1x3) : [O-O O-C C-C] temperature factors. 
% ---------------
%
% --- OUTPUTS ---
% - 'rdm' ()
% - 'dens' ()
% ---------------
% The order of the figure is the following :
% - Directional Compton Profiles (as many figures as there are directions, 3
% by default). Dashed blue curve = refined. Continuous red curve =
% pseudo-experimental.
% - Difference between previous dashed blue and continuous red ones for
% each directions.
% - Anisotropies with the last direction of the DCPs. Dashed blue curve = refined. Continuous red curve =
% pseudo-experimental. (Warning, legend has to be adapted if directions
% change, see plotCP.m file).
% - Deformation density or the difference between the refined density (diagonal terms of the 1RDM) and the pseudo-atom density map at 0K. 
% (Warning, there is no such map, like 'InputFiles/DENSPatoNew.dat', for a given temperature)
% - Refined Fourier map. Contours at +1 (blue) or -1 (red, due to truncation) x0.01x2^n (n=0 to 20) a.u^-3. 
% - Pseudo-experimental Fourier map. Contours at +1 (blue) or -1 (red, due to truncation) x0.01x2^n (n=0 to 20) a.u^-3. 
% - Difference of the previous maps : pseudo-data minus refined Fourier. Contours at 0 (green) or +1 (blue) or -1 (red) x0.0001x2^n (n=1 to 20) a.u^-3. 
% maps.
% - 1-RDM map twice along CO2 axis. Contours at +1 (blue) or -1 (red) x0.01x2^n (n=0 to 20) a.u^-3. 
% See articles for more precision on figures.
% YL.

if nargin < 9
  disp('No B factor entered');
  B =  [0 0 0];
end
%% plot CP
%plotCP(CP,P,atomicOrbitals,R);
CPpredicted = real(plotCP(CP,P,atomicOrbitals,R).');

%% map for deformation density and fourier
x=linspace(-4,4,100);
y=linspace(-3,3,75);
[map,~] = createMap(x,y);

%% plot deformation density
% Given an isolated atom map, plots the deformation density map

DENSPato = dlmread(pathPato);
plotDEFDENS(map,P,atomicOrbitals,DENSPato);

%% Plot Fourier map
% Plot Fourier map in the plane of the first CO2 molecule in the unit cell
for i=1:length(map.points)
    map.points(i,:) = (R(:,:,1)*((map.points(i,:)).')).';
end

C = orthogonalisation(atomicOrbitals);
SFOp = 0;
for i=1:size(T,1)
SFOp = SFOp + FOperators(atomicOrbitals,C,SF.Q,R(:,:,i),T(i,:));
end

SFpredicted.value = zeros(size(SFOp,3),2);
Un = ones(15,15);
for i=1:size(SFOp,3)
    DW = exp(-B*norm(SF.Q(i,:),2)^2);
    DWtemp = [DW(1)*Un DW(2)*Un Un ;
                      DW(2)*Un DW(3)*Un DW(2)*Un ;
                      Un DW(2)*Un DW(1)*Un];
    a = trace((P.*DWtemp).'*SFOp(:,:,i));
    %a = trace(P.'*SFOp(:,:,i));
    SFpredicted.value(i,1) = real(a);
    SFpredicted.value(i,2) = imag(a);
end
SFpredicted.Q = SF.Q;

plotFourier(map,SFpredicted,volume);
plotFourier(map,SF,volume);

%% Plot Fourier difference map

plotFourierError(map,SF,SFpredicted,volume);

%% Plot 1-RDM
x=-5:0.1:5;
y=-3:0.1:3;
[~,pathRDM] = createMap(x,y);

% plotRDM(pathRDM,P,atomicOrbitals);
RDMMap = plotRDM(pathRDM,P,atomicOrbitals);

set(findall(gcf,'-property','FontSize'),'FontSize',20)


end

