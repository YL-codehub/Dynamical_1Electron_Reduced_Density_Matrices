function Bfin = optimizeB(atomicOrbitals,P,R,T,SF)
%UNTITLED Summary of this function goes here
%% Orthogonalise basis set
disp(" Orthogonalizing...")
Cj = orthogonalisation(atomicOrbitals);
%% Weights
w2 = 1; % weight of Chi2_CP, to modify
%w2 = w2/(1+w2);
%w1= sqrt(1 - w2);
w2 = sqrt(w2);%weight of Chi2_SF
w1 = 1;%
%% Compute expectation operators
disp(" Computing expected operators...");
SFOp = 0;
for i=1:size(T,1)
SFOp = SFOp + FOperators(atomicOrbitals,Cj,SF.Q,R(:,:,i),T(i,:));
end

B0 = 0.001*ones(1,3); %no dilatation for O and C respectively

%%% erasing zeros
% SFzeros = ones(size(SFOp,3),1); 
% for q=1:size(SFOp,3)
%     if abs(SF.value(q,1))<=10^(-10)
%         disp('oui')
%         SFzeros(q)=0;
%     end
% end
%% fonction objectif
function val = residual(B)
    SFpred = zeros(size(SFOp,3),1);
    Un = ones(15,15);
    for q=1:size(SFOp,3)
            DW = exp(-B*norm(SF.Q(q,:),2)^2);
            DWtemp = [DW(1)*Un DW(2)*Un Un ;
                      DW(2)*Un DW(3)*Un DW(2)*Un ;
                      Un DW(2)*Un DW(1)*Un];
            SFpred(q) = trace((P.*DWtemp).'*(SFOp(:,:,q)));
    end
    %val = norm([((real(SFpred)-SF.value(:,1)).*SFzeros./SF.sigma).' ((imag(SFpred)-SF.value(:,2)).*SFzeros./SF.sigma).'], 2); %+ see later for abs for SF
  %  val = norm([((real(SFpred)-SF.value(:,1))./SF.sigma).' ((imag(SFpred)-SF.value(:,2))./SF.sigma).'], 2); %+ see later for abs for SF
   val = norm([((real(SFpred)-SF.value(:,1))./1).' ((imag(SFpred)-SF.value(:,2))./1).'], 2); %+ see later for abs for SF
end

%% declare and solve it
B = optimvar('B',length(B0),'LowerBound',0, 'Upperbound',0.03); %%%
fun = @(B) residual(B);
disp('Initial Chi^2 :');
disp(fun(B0)^2);
obj = fcn2optimexpr(fun,B);
lsqproblem = optimproblem("Objective",obj);
x0.B = B0;
show(lsqproblem);
options = optimoptions('fminunc','Display','iter-detailed', 'StepTolerance', 1e-10, 'OptimalityTolerance',1e-10);%, 'MaxFunctionEvaluations', 1000)%, 'MaxIterations',40); %see options here : https://www.mathworks.com/help/optim/ug/optimization-options-reference.html
disp(options);
[sol,fval] = solve(lsqproblem,x0,'options',options);
%% Performance :
disp('Chi2 =');
disp(fval^2);
disp('ndf = ');
disp(2*length(SF.value)-3);%2 atomes
disp('Chi2/ndf : ');
disp(fval^2/(2*length(SF.value)-3));

%% replace with new alphas

Bfin = sol.B;

%% Partial Chi^2

Chi2_SF = residual(Bfin)^2;%w1^2*norm([((real(SFpred)-SF.value(:,1))./SF.sigma).' ((imag(SFpred)-SF.value(:,2))./SF.sigma).' ],2)^2;

disp('---------------------------------------------------');
%disp('Proportion of Structure Factors in data :');
%pr = 2*length(SF.value)/(2*length(SF.value)+size(CP.value,1)*size(CP.value,2));
%disp(pr);
%disp('Proportion of Compton profiles in data :');
%disp(1-pr);

disp('---------------------------------------------------');
disp('Final Structure Factors Chi^2 :');
disp(Chi2_SF);
%disp('Final Compton profiles Chi^2 :');
%disp(Chi2_CP);

%disp('Proportion of the final Structure Factors Chi^2 :');
%disp(Chi2_SF/(fval^2));
%disp('Proportion of the final Compton profiles Chi^2 :');
%disp(Chi2_CP/(fval^2));

disp('---------------------------------------------------');
end

