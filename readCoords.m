function [ coords, R, T, volume, primVectors ] = readCoords(InName)
%% This functions aims to read the .ent file containing useful info on the crystal cell. This could also be extracted from Crystal14 outputs.
%
% --- INPUTS --- 
% - 'InName' (string) : .ent file name
% ---------------
%
% --- OUTPUTS ---
% - 'coords' (3x3): coordinates of three atoms O, C, O
% - 'R' (3x3x4): rotation matrix to rotate 'coords' to the axis of ith molecule in unit
% cell (i=1,...,4)
% - 'T' (4x3): translation vectors to go from 'coords' to the position of ith molecule in unit
% cell (i=1,...,4)
% - 'volume': volume of the primitive cell
% - 'primVectors': primitive vectors 
% ---------------
%
% YL.


fid = fopen(InName);

tline = fgetl(fid);
tline = fgetl(fid);

while ischar(tline) && startsWith(tline,'SCALE1')~=1 
      tline = fgetl(fid);
      tline;
end
SCALE = zeros(3);
for i=1:3
    line=sscanf(tline(7:end),'%f');
    SCALE(i,:)=line(1:3).';
    tline = fgetl(fid);
end
coords=zeros(2,3);
k=1;
while ischar(tline) && startsWith(tline,'CONECT')~=1 
      line=sscanf(tline(27:end-3),'%f');
      coords(k,:)=line(1:3);
      tline = fgetl(fid);
      k=k+1;
end

coords = coords*1.8897259886;

aCoords = (SCALE^-1*[1;0;0]*1.8897259886).';
bCoords = (SCALE^-1*[0;1;0]*1.8897259886).';
cCoords = (SCALE^-1*[0;0;1]*1.8897259886).';
primVectors = [aCoords;bCoords;cCoords];
volume = abs(dot(aCoords,cross(bCoords,cCoords)));


R=zeros(3,3);
T=zeros(1,3);


k=1;
for i=1:3:size(coords,1)%donc ne marche que pour molécules à 3 atomes
    av=coords(i+1,:)-coords(i,:);
    norma(k)=norm(av);

    T(k,:) = coords(i,:);

    v = cross([0;0;1],av/norma(k));
    ssc = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
    R(:,:,k) = eye(3) + ssc + ssc^2*(1-dot([0;0;1],av/norma(k)))/(norm(v))^2;
    
    k=k+1;
end
    
alength = mean(norma);

coords = [0 0 -alength; 0 0 0; 0 0 alength];
end