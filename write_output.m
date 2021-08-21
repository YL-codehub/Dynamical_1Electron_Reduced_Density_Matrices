function out = write_output(atomicOrbitals, P, T, R, OutputName)
%% This functions aims to write in an output file the resulting elements allowing to reconstruct the 1RDM. Especially, a python function exists to read it and work under python devices. 
%
% --- INPUTS --- 
% - 'atomicOrbitals' (structure) : containing as many fields as
% atomic orbitals in the model. See readBS.m for structure.
% - 'P' (matrix nOrbitals x nOrbitals) : a desorthogonalized population
% matrix. 
% - 'T' (matrix nMolecules x 3) : the list of translations allowing to go
% from reference molecule to the others in the cell.
% - 'R' (matrix nMolecules x 3x3) : the list of rotation matrices allowing to go 
% from reference molecule to the others in the cell.
% - 'OutputName' (string) : the name of the file we write in.
% ---------------
%
% --- OUTPUTS ---
% - 'out' : /
% ---------------
%
% YL.

disp('Please pay attention to desorthogonalize your matrix. If solP is your matrix coming from SDP algorithm, and C your orthonalizing matrix, you have to write: C^T x solP x C');

%% Write atomic orbitals in specific format
file = fopen(OutputName, 'w');
fprintf(file, 'LCAO');
fprintf(file, '\n');
for i = 1:size(atomicOrbitals,2)
    temp = atomicOrbitals{1, i}.d;
    for k = 1:size(temp,2)
        fprintf(file,string(temp(1,k)));
        fprintf(file, '\t');
    end
    fprintf(file, ' |');
    for j = 1:size(atomicOrbitals{1, i}.primitives,2)
        gto = atomicOrbitals{1, i}.primitives{1, j};
        fprintf(file, '; ');
        fprintf(file,string(gto.a(1))); fprintf(file, '\t'); fprintf(file,string(gto.a(2))); fprintf(file, '\t'); fprintf(file,string(gto.a(3))) ; %a
        fprintf(file, ' & ');
        fprintf(file,string(gto.A(1))); fprintf(file, '\t'); fprintf(file,string(gto.A(2)));fprintf(file, '\t'); fprintf(file,string(gto.A(3))) ; %A
        fprintf(file, ' & ');    
        fprintf(file,string(gto.alpha));
    end
    fprintf(file, '\n');
end

fprintf(file, 'Population matrix');
fclose(file);

%% Write population matrix
writematrix(P, OutputName,'Delimiter', 'tab','WriteMode','append');

%% Write translations and rotations
file = fopen(OutputName, 'a+');
fprintf(file, 'Translations');
fclose(file);
writematrix(T, OutputName,'Delimiter', 'tab','WriteMode','append');

file = fopen(OutputName, 'a+');
fprintf(file, 'Rotations');
fclose(file);
for i=1:size(R,3)
    writematrix(R(:,:,i), OutputName,'Delimiter', 'tab','WriteMode','append');
end