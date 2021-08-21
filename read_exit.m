function [atomicOrbitals, P, R, T] = read_exit(FileName,gamess)
if nargin < 2
  gamess = false;
end
% normalement l'ordre des orbitales est conservée
page = fopen(FileName);
tline = fgetl(page);
tline = fgetl(page);
atomicOrbitals = {};
while ischar(tline) && startsWith(tline,'Population matrix')~=1 %=before pop matrix
    a = split(tline, '|;');
    contractionCoefs=sscanf(a{1}, '%f').';
    aos = split(a{2}, ';');
    orbType = ' ';
    powerExponents=[];
    atomCount = 1; %useless..
    atom = 0; %lack of atomic numbers in exitmatlab files but useless
    nShell = 1;
    coords = [];
    for i=1:length(aos)
        caracs = split(aos{i}, '&');
        L = sum(sscanf(caracs{1}, '%f'));
        powerExponents = [powerExponents sscanf(caracs{3}, '%f')]; % /!\ %f arrondit les valeurs...
    end
    coords =  sscanf(caracs{2},'%f').'; %not necessary to check for each gaussian because it's the same for all
    if L == 0 %not necessary to check for each gaussian because it's the same for all
        orbType = 's';
    elseif L == 1
        orbType = 'p';
    else
        orbType = 'd';
    end
    %atomicShells{atomCount,nShell}=atomicShell(coords,atom,orbType,contractionCoefs,powerExponents);
    atomicOrbitals = {atomicOrbitals{:}, atomicOrbital(coords,0,sscanf(caracs{1}, '%f').',contractionCoefs,powerExponents)};
    tline = fgetl(page);
    nShell = nShell+1;
end
tline = fgetl(page);
P = [];
R_ = [];
T = [];
 while ischar(tline) && startsWith(tline,'Translations')~=1 %before translations
     P = [ P ; sscanf(tline, '%f').'];
     tline = fgetl(page);
 end

%  %% Gamess' case : Reverse orbitales d
function out = permute(a,i,j)
            b = a(:,j);
            a(:,j)=a(:,i);
            a(:,i)=b;
            b = a(j,:);
            a(j,:)=a(i,:);
            a(i,:)=b;
            out = a;
end
     
if gamess
    P = permute(P,11,13);
    P = permute(P,12,14);
    P = permute(P,14,15);
    P = permute(P,26,28);
    P = permute(P,27,29);
    P = permute(P,29,30);
    P = permute(P,41,43);
    P = permute(P,42,44);
    P = permute(P,44,45);

    for i=1:15
        atomicOrbitals{1, i}.atomType = 6;  
        atomicOrbitals{1, i+15}.atomType = 8;
        atomicOrbitals{1, i+30}.atomType = 8;
        P = permute(P,i,i+15);
        temp = atomicOrbitals{1,i};
        atomicOrbitals{1,i}=atomicOrbitals{1,i+15};
        atomicOrbitals{1,i+15}=temp;
    end
end

C = orthogonalisation(atomicOrbitals);
P = C.'^-1*P*C^-1;
 
while ischar(tline) && startsWith(tline,'Rotations')~=1 %before Rotations
     T = [ T ; sscanf(tline, '%f').'];
     tline = fgetl(page);
end

while ischar(tline)  %before end of file
     R_ = [ R_ ; sscanf(tline, '%f').'];
     tline = fgetl(page);
end

R = zeros(3,3,length(R_)/3);
for i=1:length(R_)/3
   R(:,:,i)= R_(1+3*(i-1):3+3*(i-1),:);
end
end