function [SF,CP] = readCrystal(primVectors,InputPropName,FSname,CPname, Rinf,Rsup, qM,step)
%% This functions aims to read the data points from Crystal 14 Software outputs.
%
% --- INPUTS --- 
% - 'primVectors' (string) : primitive vectors of the crystal cell given by 'readCoords' function.
% - 'InputPropName' (string) : name of the input file for the computation of SFs and
% DCPs (properties input file, usually with .d3 extension) containing data directions asked to Crystal.
% - 'FSname' (string) : output file of properties computation by Crystal (usually
% with .out extension). It is used to read SFs.
% - 'CPname' (string) :  output file of DCPs computation by Crystal (usually with a
% .DAT extension and named 'PROF')
% - 'Rinf' (float) : minimal radius of the Ewald sphere of data which is read for
% SFs. For dry ice, a radius 8 gives a sin(theta)/lambda = 0.71.
% - 'Rsup' (float) :  maximal radius of the Ewald sphere of data which is read for
% SFs. For dry ice, a radius 12 gives a sin(theta)/lambda = 1.07.
% - 'qM' (float) : maximal value of momentum which is read for DCPs.
% - 'step' (float) : the DCPs resolution or the impulsion difference between two
% points of a DCP. It is usually no less than 0.1 a.u.
% ---------------
%
% --- OUTPUTS ---
% - 'SF' (structure) : containing SFs thanks to 3 fields :
%            - SF.Q (nSFs x 3): the list of directions in Angströms. The
%            number nSFs of directions, depends on Rinf and Rsup and of the
%            work of Crystal since it usually removes equivalent
%            directions.
%            - SF.value (nSFs x 2) : the list of the associated values, 1st
%            column being the real part and the second the imaginary one.
%            The first value (corresponding to [0 0 0]) gives the number of electrons in the cell. 
%            - SF.sigma (nSFs x 1): the associated errorbars. Set to 1 by
%            default since pseudo-data from Crystal are noise-free.
% - 'CP' (structure) : containing DCPs thanks to 4 fields :
%            - CP.e (nDCPs x 3): the list of directions in Angströms. The
%            number nDCPs of directions is chosen before Crystal calculus
%            in the input file.
%            - CP.q (nSteps x 1) : the list of impulsion points considered
%            for all directions. It is given by the choice of resolution
%            since nSteps = qM/step.
%            - CP.value (nSteps x nDCPs) : the list of the associated
%            values, each column being the value for each direction.
%            - CP.sigma (nSteps x nDCPs): the associated errorbars. Set to 1 by
%            default since pseudo-data from Crystal are noise-free.
% ---------------
%
% YL.

%% Reading Structure factors 

fid = fopen(FSname);
tline = fgetl(fid);

while ischar(tline) && startsWith(tline,' NUMBER OF FACTORS COMPUTED :')~=1 
      tline = fgetl(fid);
end
tline = regexprep(tline,':',' '); tline = regexprep(tline,'\w*\s',' ');
Nqq = sscanf(tline,'%d').';

while ischar(tline) && startsWith(tline,'    ALPHA+BETA ELECTRONS')~=1 
      tline = fgetl(fid);
end

fgetl(fid);fgetl(fid);tline=fgetl(fid);
 
k = 1 ;%
hkl = zeros(Nqq,3);
SFexp = zeros(Nqq,2);

%data selection in Ewald sphere between radius Rinf and Rsup  
for i=1:Nqq
      line = sscanf(tline,'%f').';
      r = (line(1)^2 + line(2)^2 + line(3)^2);
      if r<=Rsup^2 && r>=Rinf^2 %& SF_hole(line(1:3),Rsup)%uncomment to hole data
        hkl(k,:) = line(1:3);
        SFexp(k,:) = line(end-2:end-1);
        k = k+1;
      end
      tline = fgetl(fid);
end


hkl = hkl(1:k-1,:);
SFexp(1:k-1,:);
fclose(fid);

%% Uncomment and use this code if you want to read a parallepideic data set as B. De Bruyne did
% hkl = zeros(Nqq,3);
% SFexp = zeros(Nqq,2);
% for i=1:Nqq
%       line = sscanf(tline,'%f').';
%       hkl(i,:) = line(1:3);
%       SFexp(i,:) = line(end-2:end-1);
%       tline = fgetl(fid);
% end
%lim = Rsup; %11.3;%Benjamin took lim=7
% for i=1:Nqq
%       line = sscanf(tline,'%f').';
%       if 0<=line(1) & line(1)<=lim & -lim<=line(2) & line(2)<=lim & -lim<=line(3) & line(3)<=lim %
%         hkl(k,:) = line(1:3);
%         SFexp(k,:) = line(end-2:end-1);
%         k = k+1;
%       end
%       tline = fgetl(fid);
% end


%% Reading Directionnal Compton Profiles

% Read directions
fid = fopen(InputPropName);
tline = fgetl(fid);

while ischar(tline) && startsWith(tline,'CP')~=1
      tline = fgetl(fid);
end
tline=fgetl(fid);
line = sscanf(tline,'%f').';
Nhkl = line(1);
qMax = line(2);

CPehkl = zeros(Nhkl,3);
for i=1:Nhkl
    tline = fgetl(fid);
    CPehkl(i,:)=sscanf(tline,'%f').';
end

tline=fgetl(fid);
line = sscanf(tline,'%f');
CPvaluesQ=0:line(1):qMax;
Ncp=length(CPvaluesQ);

fclose(fid);

% read values
fid = fopen(CPname);
tline = fgetl(fid);

CPexp=zeros(Nhkl,Ncp);
for i=1:Nhkl
    for j=1:Ncp
        line = sscanf(tline,'%f').';
        CPexp(i,j) = line(2);
        tline = fgetl(fid);
    end
    fgetl(fid);tline=fgetl(fid);
end

fclose(fid);

%% Convert vectors from Ewald sphere to real space
a1 = primVectors(1,:);
a2 = primVectors(2,:);
a3 = primVectors(3,:);

b1 = 2*pi*cross(a2,a3)/(dot(a1,cross(a2,a3)));
b2 = 2*pi*cross(a3,a1)/(dot(a2,cross(a3,a1)));
b3 = 2*pi*cross(a1,a2)/(dot(a3,cross(a1,a2)));

Q = zeros(length(hkl),3);
for i=1:length(hkl)
    Q(i,:) = hkl(i,1)*b1 + hkl(i,2)*b2 + hkl(i,3)*b3; 
end

e = zeros(length(CPehkl),3);
for i=1:length(CPehkl)
    e(i,:) = CPehkl(i,1)*b1 + CPehkl(i,2)*b2 + CPehkl(i,3)*b3; 
    e(i,:) = e(i,:)/norm(e(i,:));
end


%% Constrain and store values
qMax= qM;
pas = step;
indQMax = find(CPvaluesQ <=qMax,1, 'last');
CPexp=CPexp(:,1:pas:indQMax);
CPvaluesQ=CPvaluesQ(1:pas:indQMax);

Nqq=size(Q,1);
Ncp=length(CPvaluesQ);
SF.Q = Q(1:Nqq,:);
SF.value = SFexp(1:Nqq,:);
CP.e = e;
CP.q = CPvaluesQ(1:Ncp).';
CP.value = CPexp(:,1:Ncp).';

%% Set default error bars
sigmaSF = (SF.value(:,1)).^0;
sigmaCP = sqrt(CP.value).^0;

SF.sigma = sigmaSF;
CP.sigma = sigmaCP;

end