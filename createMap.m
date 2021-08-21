function [map,RDMpath] = createMap(x,y)


[X,Y] = meshgrid(x,y);
map.points = zeros(length(x)*length(y),3);
map.N = length(y);
k=1;
for j=1:length(x)
    for i=1:length(y)
        map.points(k,:) = [0,Y(i,j),X(i,j)];
        k=k+1;
    end
end

x=-3:0.1:3;
RDMpath = [0*ones(length(x),1),zeros(length(x),1),x.'];
%RDMpath = [x.',0*ones(length(x),1),zeros(length(x),1)];
end