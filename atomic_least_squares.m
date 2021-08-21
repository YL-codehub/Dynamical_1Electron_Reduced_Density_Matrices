function atomicOrbitals_final = atomic_least_squares(atomicOrbitals,P,R,T,SF,CP,B)
%% This function aims to contract or dilate the basis set with a least squares method on entered data, with or without deconvolution of SFs thanks to temperature factor B.
%
% --- INPUTS --- 
% - 'P' (nOrbs x nOrbs) : Population matrix from refinements (still
% orthogonalized)
% - 'atomicOrbitals' (structure) : initial non-orthogonalized atomic orbitals.
% - 'R' (nMolecules x 3x3) : rotation matrices to go from reference molecules
% to all in the cell.
% - 'T' (nMolecules x 3) : DCPs data points.
% - 'SF' (structure) : SFs data points.
% - 'CP' (structure) : DCPs data points.
% - 'B' (1x3) : [O-O O-C C-C] temperature factors. [0 0 0] if not entered.
% ---------------
%
% --- OUTPUTS ---
%- 'atomicOrbitals' (structure) : contracted-dilated non-orthogonalized atomic orbitals.
% ---------------
%  The same coefficient is applied for every GTO of the same atom. To
%  contract/dilate only valence GTOs, see 'valence_least_squares.m' in
%  'Explorations' file.
% YL.

if nargin < 7
  disp('No B factor entered');
  B =  [0 0 0];
end

%% Set Chi^2 Weights
w2 = 1; % weight of Chi2_CP, to be modified
%w2 = w2/(1+w2);
%w1= sqrt(1 - w2);
%w2 = sqrt(w2);%weight of Chi2_SF
w1 = 1;%

%% Initial solver guess point
kappas0 = [0.95 0.95]; %no dilatation for O and C respectively

%% Residual function
function val = residual(Kappas)
%disp(Kappas);
atomicOrbitals_new = atomicOrbitals;
for i = 1:length(atomicOrbitals_new)
    if i<=15 || i>= 31
        for j=1:length(atomicOrbitals_new{1, i}.d)
           atomicOrbitals_new{1, i}.primitives{1, j}.alpha = Kappas(1)*atomicOrbitals_new{1, i}.primitives{1, j}.alpha; % à modifier
        end
    else
        for j=1:length(atomicOrbitals_new{1, i}.d)
           atomicOrbitals_new{1, i}.primitives{1, j}.alpha = Kappas(2)*atomicOrbitals_new{1, i}.primitives{1, j}.alpha; % à modifier
        end
    end
               
end
[SFpred, CPpred] = predictor(P,atomicOrbitals_new, R, T, SF, CP,B);
val = norm([((real(SFpred)-SF.value(:,1)).*w1./SF.sigma).' ((imag(SFpred)-SF.value(:,2)).*w1./SF.sigma).' reshape((CPpred.'-CP.value).*w2./CP.sigma,[1,size(CP.value,1)*size(CP.value,2)])  ], 2); %+ see later for abs for SF
end
%% declare and solve it
Kappas = optimvar('Kappas',length(kappas0),'LowerBound',0,'UpperBound',2); %%%
fun = @(Kappas) residual(Kappas);
disp('Initial Chi^2 :');
disp(fun(kappas0)^2);
obj = fcn2optimexpr(fun,Kappas);
lsqproblem = optimproblem("Objective",obj);
x0.Kappas = kappas0;
show(lsqproblem);
options = optimoptions('fminunc','Display','iter-detailed');%, 'StepTolerance', 1e-10, 'OptimalityTolerance',1e-10)%, 'MaxFunctionEvaluations', 1000)%, 'MaxIterations',40); %see options here : https://www.mathworks.com/help/optim/ug/optimization-options-reference.html
disp(options);
[sol,fval] = solve(lsqproblem,x0,'options',options);

%% Performance :
disp('Chi2 =');
disp(fval^2);
disp('ndf = ');
disp(2*length(SF.value)+size(CP.value,1)*size(CP.value,2)-2);%2 atomes
disp('Chi2/ndf : ');
disp(fval^2/(2*length(SF.value)+size(CP.value,1)*size(CP.value,2)-2));

%% Replace with new alphas

kappas = sol.Kappas;
disp(kappas);
atomicOrbitals_final = atomicOrbitals;
for i = 1:length(atomicOrbitals_final)
    if i<=15 || i>= 31
        for j=1:length(atomicOrbitals_final{1, i}.d)
           atomicOrbitals_final{1, i}.primitives{1, j}.alpha = kappas(1)*atomicOrbitals_final{1, i}.primitives{1, j}.alpha; % à modifier
        end
    else
        for j=1:length(atomicOrbitals_final{1, i}.d)
           atomicOrbitals_final{1, i}.primitives{1, j}.alpha = kappas(2)*atomicOrbitals_final{1, i}.primitives{1, j}.alpha; % à modifier
        end
    end
               
end

%% Partial Chi^2
[SFpred, CPpred] = predictor(P,atomicOrbitals_final, R, T, SF, CP);
% Chi2_SF = norm([((real(SFpred)-SF.value(:,1))./SF.sigma).' ((imag(SFpred)-SF.value(:,2))./SF.sigma).' ],2)^2;
% Chi2_CP = norm([reshape((CPpred.'-CP.value)./CP.sigma,[1,size(CP.value,1)*size(CP.value,2)]) ], 2)^2;

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
disp('Final Compton profiles Chi^2 :');
disp(Chi2_CP);

disp('Proportion of the final Structure Factors Chi^2 :');
disp(Chi2_SF/(fval^2));
disp('Proportion of the final Compton profiles Chi^2 :');
disp(Chi2_CP/(fval^2));

disp('---------------------------------------------------');
end