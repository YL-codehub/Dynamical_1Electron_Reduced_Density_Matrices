function [solP, Cj] = optimizeP_quadratic(atomicOrbitals,R,T,SF,CP,B)
%% This function aims to refine the population matrix on entered data, with or without deconvolution of SFs thanks to temperature factor B.
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
% - 'B' (1x3) : [O-O O-C C-C] temperature factors. [0 0 0] if not entered.
% ---------------
%
% --- OUTPUTS ---
% - 'solP' (nOrbs x nOrbs) : Population matrix from refinements (still
% orthogonalized)
% - 'Cj' (nOrbs x nOrbs) : Lowdin basis change matrix for orthog/
% desorthogonalization.
% ---------------
% 
% YL.

if nargin < 6
  disp('No B factor entered');
  B =  [0 0 0];
end

%% P as unknown
P = sdpvar(length(atomicOrbitals));

%% Orthogonalise basis set
disp(" Orthogonalizing...")
Cj = orthogonalisation(atomicOrbitals);

%% Compute expectation operators
disp(" Computing expected operators...");
SFOp = 0;
CPOp = 0;
for i=1:size(T,1)
SFOp = SFOp + FOperators(atomicOrbitals,Cj,SF.Q,R(:,:,i),T(i,:));
CPOp = CPOp + CPOperators(atomicOrbitals,Cj,CP.e,CP.q,R(:,:,i));
end

%% Compute expectation values
disp(" Creating Mosek problem instance :");
disp(" Declaring SF predictor ");
SFpredicted = sdpvar(size(SFOp,3),1);
Un = ones(15,15);
Zer = zeros(15,15);

for i=1:size(SFOp,3)
    DW = exp(-B*norm(SF.Q(i,:),2)^2);
    DWtemp = [DW(1)*Un DW(2)*Un Un ;
               DW(2)*Un DW(3)*Un DW(2)*Un ;
                  Un DW(2)*Un DW(1)*Un];
    SFpredicted(i) = trace((P.*DWtemp).'*SFOp(:,:,i));
end
disp(" Declaring CP predictor ");

CPpredicted = sdpvar(size(CPOp,3)*size(CPOp,4),1);
k=1;
for e=1:size(CPOp,3)
    for q=1:size(CPOp,4)
        CPpredicted(k) = trace(P.'*CPOp(:,:,e,q));
        k=k+1;
    end
end

%% Difference between computed and given
disp(" Chi2 declaration...");

%w2 = 1; % weight of Chi2_CP, to modify
%w2 = w2/(1+w2);
%w1= sqrt(1 - w2);


wit = 1;
for wi=wit:1:wit
    w2 = sqrt(wi);%weight of Chi2_SF
    w1 = 1;%
   
    Difference = [...
        (real(SFpredicted)-SF.value(:,1)).*w1./SF.sigma;...
        (imag(SFpredicted)-SF.value(:,2)).*w1;...
        (real(CPpredicted)-CP.value(:)).*w2./CP.sigma(:)...
    ];
  
    %% Constraints
    disp(" Declaring constraints...");
    t = norm(Difference,2);
    
    C = [ P>= 0];
    C = [C, 2*eye(length(P)) - P >= 0];
    C = [C, trace(P)== SF.value(1,1)/size(T,1)];

    %% Optimization
    disp(" Optimization in progress...");
    sol = optimize(C,t)
    %% Partial Chi^2
    [SFpred, CPpred] = predictor(value(P),atomicOrbitals, R, T, SF, CP,B);
    
    Chi2_SF = w1^2*norm([((real(SFpred)-SF.value(:,1))./SF.sigma).' ((imag(SFpred)-SF.value(:,2))./SF.sigma).' ],2)^2;
    Chi2_CP = w2^2*norm([reshape((CPpred.'-CP.value)./CP.sigma,[1,size(CP.value,1)*size(CP.value,2)]) ], 2)^2;

    disp('---------------------------------------------------');
    disp('Proportion of Structure Factors in data :');
    pr = 2*length(SF.value)/(2*length(SF.value)+size(CP.value,1)*size(CP.value,2));
    disp(pr);
    disp('Proportion of Compton profiles in data :');
    disp(1-pr);

    disp('---------------------------------------------------');
    disp('Final Structure Factors Chi^2 :');
    disp(Chi2_SF);
    disp('Final weighted Compton profiles Chi^2 :');
    disp(Chi2_CP);
    disp('Final Compton profiles Chi^2 :');
    disp(Chi2_CP/w2^2);

    disp('Proportion of the final Structure Factors Chi^2 :');
    disp(Chi2_SF/(value(t)^2));
    disp('Proportion of the final weighted Compton profiles Chi^2 :');
    disp(Chi2_CP/(value(t)^2));
    disp('Proportion of the final Compton profiles Chi^2 :');
    disp(Chi2_CP/w2^2/(Chi2_SF+Chi2_CP/w2^2));

    disp('---------------------------------------------------');
    %% Quality of optimization
    disp('Chi2 =');
    disp(value(t)^2);
    disp('ndf = ');
    disp(2*length(SF.value)+size(CP.value,1)*size(CP.value,2)-(0.5*size(P,1)*(size(P,2)+1)-1));
    disp(wi);
    disp('Chi2/ndf : ');
    disp(value(t)^2/(2*length(SF.value)+size(CP.value,1)*size(CP.value,2)-(0.5*size(P,1)*size(P,2)-1)));
end
%% Solution
solP = value(P);
end

