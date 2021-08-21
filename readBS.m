function [atomicOrbitals,M] = readBS(InputName,coords)
%% This functions aims to read the modeling basis set.
%
% --- INPUTS --- 
% - 'InputName' (string) : name of the .d12 input file (Crystal format)
% containing the basis set which will be used for the refined elements (model).
% - 'coords' (matrix 3 x 3) : coordinates of atoms in the reference CO2
% molecule (those are rotated and translated later to get the 4 molecules of
% a cell in the different refinements).
% ---------------
%
% --- OUTPUTS ---
% - 'atomicOrbitals' (structure) : containing as many fields as
% atomic orbitals in the model. Note that f orbitals cannot be read. One atomic orbital i is called with
% atomicOrbitals{1,i} and contains (see atomicOrbital class file):
%                 - atomicOrbitals{1,i}.d : a line with the contraction
%                 coefficients of the primitive functions (see below).
%                 - atomicOrbitals{1,i}.primitives : a structure containing
%                 as much fields as there are gaussian primitive functions
%                 in this atomic orbital. Such a primitive contains the fields : 
%                              - atomicOrbitals{1, i}.primitives{1, j}.A : coordinates of the gaussian function
%                              - atomicOrbitals{1, i}.primitives{1, j}.a : exponents of (X-A_X), X = x,y or z
%                              - atomicOrbitals{1, i}.primitives{1,j}.alpha : exponent of the gaussian
%                              function.
%                              - atomicOrbitals{1, i}.primitives{1, j}.N : normalization coefficient 
%                 - atomicOrbitals{1,i}.atomType : the atomic number of the
%                 atom to which the orbital is attached.
% - 'M' : passage matrix from spherical to cartesian coordinates. In fact d
% or f atomic orbitals are often purely spherical in Crystal software. This
% is why a basis chage is necessary in this case. (see atomicShell class file)
%
% YL.


fid = fopen(InputName);
tline = fgetl(fid);


while ischar(tline) && startsWith(tline,'ENDG')~=1 
      tline = fgetl(fid);
end
tline = fgetl(fid);

atomCount=1;
while ischar(tline) && startsWith(tline,'99')~=1 
    
    line = sscanf(tline,'%f');
    atom = line(1);
    numShell = line(2);
    for nShell=1:numShell
        powerExponents=[];
        contractionCoefs=[];
        tline = fgetl(fid);
        line = sscanf(tline,'%f');
        type = line(2);
        numG = line(3);

        if type==0
            orbType = 's';
            for g=1:numG
                tline = fgetl(fid);
                line = sscanf(tline,'%f');
                powerExponents(g)=line(1);
                contractionCoefs(g)=line(2);
            end
            atomicShells{atomCount,nShell}=atomicShell(coords(atomCount,:),atom,orbType,contractionCoefs,powerExponents);
            
         elseif type==1
             orbType = 'sp';
            for g=1:numG
                tline = fgetl(fid);
                line = sscanf(tline,'%f');
                powerExponents(g)=line(1);
                contractionCoefs(1,g)=line(2);
                contractionCoefs(2,g)=line(3);
            end
            atomicShells{atomCount,nShell}=atomicShell(coords(atomCount,:),atom,orbType,contractionCoefs,powerExponents);
            
        elseif type==2
            orbType = 'p';
            for g=1:numG
                tline = fgetl(fid);
                line = sscanf(tline,'%f');
                powerExponents(g)=line(1);
                contractionCoefs(g)=line(2);
            end
            atomicShells{atomCount,nShell}=atomicShell(coords(atomCount,:),atom,orbType,contractionCoefs,powerExponents);
         elseif type==3
            orbType = 'd';
            for g=1:numG
                tline = fgetl(fid);
                line = sscanf(tline,'%f');
                powerExponents(g)=line(1);
                contractionCoefs(g)=line(2);
            end
            atomicShells{atomCount,nShell}=atomicShell(coords(atomCount,:),atom,orbType,contractionCoefs,powerExponents);
        end
    end
    
    tline=fgetl(fid);
    atomCount=atomCount+1;
end

%% Store atomic orbitals thanks to atomicShells and atomic atomicOrbitals classes in Cartesian coordinates
atomicOrbitals={};
M = [];
for atomCount=1:size(atomicShells,1)
    for shell=1:size(atomicShells,2)
        if ~isempty(atomicShells{atomCount,shell})
        atomicOrbitals={atomicOrbitals{:},atomicShells{atomCount,shell}.atomicOrbitals{:}};
        M = blkdiag(M,atomicShells{atomCount,shell}.CartesianToSpherical);
        end
    end
end

end