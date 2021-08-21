function CPpredicted = plotCP(CP,P,atomicOrbitals,R)
% function gca = plotCP(CP,P,atomicOrbitals,R)
%% This functions aims to plot DCPs performances : refined vs data ones, differences and anisotropies.
%
% --- INPUTS --- 
% - 'CP' (structure) : DCPs data points (see readCrystal).
% - 'P' (nOrbs x nOrbs) : Population matrix from refinements (still
% orthogonalized)
% - 'atomicOrbitals' (structure) : non-orthogonalized atomic orbitals.
% - 'R' (nMolecules x 3) : rotation matrices to go from reference molecules
% to all in the cell.
% ---------------
%
% --- OUTPUTS ---
% - 'gca' (plot)
% ---------------
%
% YL.

C = orthogonalisation(atomicOrbitals);

%% Compute expectation operators
CPOp = 0;
for i=1:size(R,3)
CPOp = CPOp + CPOperators(atomicOrbitals,C,CP.e,CP.q,R(:,:,i));
end

%% Compute expectation values
CPpredicted = zeros(size(CPOp,3),size(CPOp,4));
for e=1:size(CPOp,3)
    for q=1:size(CPOp,4)
        CPpredicted(e,q) = trace(P.'*CPOp(:,:,e,q));
    end
end

R = CP;
R.value = sum(abs(abs(R.value)-abs(CPpredicted.')),[1 2])/sum(abs(R.value),[1 2]);
disp('R factor for DCPs is :');
R.value

for i=1:size(CP.e,1)
    figure;
%     gca = plot(CP.q,(CP.value(:,i))/(trace(P)*size(R,3)),'r',CP.q,real(CPpredicted(i,:))/(trace(P)*size(R,3)),'b','LineWidth',1);
    plot(CP.q,(CP.value(:,i))/(trace(P)*size(R,3)),'r','LineStyle','-','LineWidth',1);
    hold on;
    plot(CP.q,real(CPpredicted(i,:))/(trace(P)*size(R,3)),'b','LineStyle','--','LineWidth',1);
    ylabel('J(q) (a.u.)');
  % legend('Crystal','Predicted');
    title(['R = ', num2str(R.value)],'interpreter','latex')
    xlabel(['q (a.u.)']);
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
end
% J'ai l'impression que la normalisation est cassée... voir tracé qui n'est
% plus autour de 0.3

%% Difference CP - CPexp

for i=1:size(CP.e,1)
    figure;
    temp = (CP.value(:,i).'-real(CPpredicted(i,:)));
    %temp = (CP.value(:,i).'-real(CPpredicted(i,:)))./CP.value(:,i).';
    gca = plot(CP.q,temp/(trace(P)*size(R,3)),'r','LineWidth',2);
    %gca = plot(CP.q,temp,'r','LineWidth',2);
    ylabel('(Jexp(q)-J(q)) (a.u.)');
  % legend('Crystal','Predicted');
    xlabel(['q (a.u.)']);
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
end
set(findall(gcf,'-property','FontSize'),'FontSize',20);

%% Anisotropies
figure;
labels = [0 0 1 ; 1 1 0 ; 1 1 1];

for i=1:size(CP.e,1)-1
    subplot(1,size(CP.e,1)-1,i);
    set(gcf, 'Position',  [500, 500, 900, 300]);
    plot(CP.q,CP.value(:,i)/(trace(P)*size(R,3))-CP.value(:,3)/(trace(P)*size(R,3)),'r','LineStyle','-','LineWidth',1);
    hold on;
    plot(CP.q,real(CPpredicted(i,:))/(trace(P)*size(R,3))-real(CPpredicted(3,:))/(trace(P)*size(R,3)),'b','LineStyle','--','LineWidth',1);
    title(['[ ', num2str(labels(i,:)),' ]-[', num2str(labels(3,:)),' ]'],'interpreter','latex')
  %  legend('Crystal','Predicted');
    xlabel(['q( a.u. )'],'interpreter','latex')
    ylabel(['$$\Delta J(q)$$ (a.u.)'],'interpreter','latex');
end
sgtitle(['$$R_J = ', num2str(R.value), '$$'],'interpreter','latex');


end

