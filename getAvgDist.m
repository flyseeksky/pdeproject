function [avgDist] = getAvgDist(points, n)
sum = 0;
for i = 1:n
    sum = sum + norm(points(1:2, i) - points(1:2, mod(i,n)+1));
end
avgDist = sum/n;
end
