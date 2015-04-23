function [newPoints] = snakeResample(points)

% Get number of points in the snake
n = size(points, 2);

% Get number of rows in the points matrix
rows = size(points, 1);

% Get the average distance between snake points
avgDist = getAvgDist(points, n);

% Calculate distance between each pair of points
dist = zeros(n,1);
for i = 1:n
    dist(i) =  norm(points(1:2, i) - points(1:2, mod(i,n)+1));
end

% Find places in the snake where additional control points
% needs to be inserted
insertAt = find(dist >= 1.3*avgDist);

% Number of new points that should be added
add = length(insertAt);

% Find places in the snake where control points
% needs to be removed
removeAt = find(dist <= 0.6*avgDist);

% Number of points that should be removed
remove = length(removeAt);

% Calculate number of points in the resampled snake
newPoints = zeros(rows, n+add-remove);
m = 0;
i = 1;

% Loop that does the actual resampling
while i <= n
    [modI, modIminus, modIplus] = getModulo(i, n);
    if any(insertAt == i)
        newPoints(:, i+m) = points(:, i);
        newPoints(1, i+1+m) = ( points(1, modIplus) - points(1, modI) )/2 + points(1,modI);
        newPoints(2, i+1+m) = ( points(2, modIplus) - points(2, modI) )/2 + points(2,modI);
        if rows == 3
            newPoints(3, i+1+m) = ( points(3, modIplus) + points(3, modI) )/2;
        elseif rows == 4
            newPoints(3, i+1+m) = ( points(3, modIplus) + points(3, modI) )/2;
            newPoints(4, i+1+m) = ( points(4, modIplus) + points(4, modI) )/2;
        elseif rows == 5
            newPoints(3, i+1+m) = ( points(3, modIplus) + points(3, modI) )/2;
            newPoints(4, i+1+m) = ( points(4, modIplus) + points(4, modI) )/2;
            newPoints(5, i+1+m) = ( points(5, modIplus) + points(5, modI) )/2;
        end
        m = m + 1; 
    else
        newPoints(:, i+m) = points(:, i);
    end
    % Remove points
    if any(removeAt == i)
        newPoints(:, i+1+m) = points(:, modI);
        m = m - 1; 
    end
    i = i + 1; 
end

end

