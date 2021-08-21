function [SFpred, CPpred] = predictor(P,atomicOrbitals, R, T, SF, CP, B)
%% This functions aims to provide the list of values of the modelized SFs and DCPs data points. 
% This reproduces the calculus made inside optimizers but independantly in order to compute Chi2/ndf performances. 
%
% --- INPUTS --- 
% - 'P' (nOrbs x nOrbs) : Population matrix from refinements (still
% orthogonalized)
% - 'atomicOrbitals' (structure) : non-orthogonalized atomic orbitals.
% - 'R' (nMolecules x 3) : rotation matrices to go from reference molecules
% to all in the cell.
% - 'SF' (structure) : SFs data points (see readCrystal).
% - 'CP' (structure) : DCPs data points (see readCrystal).
% - 'B' (1x3) : [O-O O-C C-C] temperature factors. 
% ---------------
%
% --- OUTPUTS ---
% - SFpred (nSF.Q x 2) : the list of real (1st column) and imaginary (2nd
% column) parts of the modelized structure factors. Each line corresponds
% to the Q-vector in the same line in SF.Q
% - CPpred (nCP.e x nCP.q     x 1) :  the list of values of the modelized
% DCPs. Each third of the list corresponds to a direction in CP.e, 
% each line in such a third corresponds to the line in CP.Q.
% ---------------
%
% YL.


%% if no temperature factor:
if nargin < 7
  disp('No B factor entered');
  B =  [0 0 0];
end

%% Orthogonalization is necessary
Cj = orthogonalisation(atomicOrbitals);

%% Compute operators
SFOp = 0;
CPOp = 0;
for i=1:size(T,1)
SFOp = SFOp + FOperators(atomicOrbitals,Cj,SF.Q,R(:,:,i),T(i,:));
CPOp = CPOp + CPOperators(atomicOrbitals,Cj,CP.e,CP.q,R(:,:,i));
end

%% Compute expectation values
SFpred = zeros(size(SFOp,3),1);
for q=1:size(SFOp,3)
        DW = exp(-B*SF.Q(q)^2);
        Un = ones(15,15);
        DWtemp = [DW(1)*Un DW(2)*Un Un ;
                      DW(2)*Un DW(3)*Un DW(2)*Un ;
                      Un DW(2)*Un DW(1)*Un];
        SFpred(q) = trace((P.*DWtemp).'*SFOp(:,:,q));
end

CPpred= zeros(size(CPOp,3),size(CPOp,4));
for e=1:size(CPOp,3)
    for q=1:size(CPOp,4)
        CPpred(e,q) = trace(P.'*CPOp(:,:,e,q));
    end
end

end
