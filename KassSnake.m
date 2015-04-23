function [cellArray] = KassSnake(points, img, alpha, beta, delta, sigma, maxIt, rs)

    x = points(1,:)';
    y = points(2,:)';

% Thresholds used in terminating the algorithm
stopThresh = 0.0008;
stop = 0;
scale = 1;

% Get number of points in the snake
n = size(points, 2);

% Preallocations for increased code efficiency
cellArray = cell(maxIt, 1);

% Get the image energy
enrgImg = -getImgEnrg(img, sigma);

% Normalize image energy by Min-Max normalization
% A = (A - minA)/(maxA - minA) * (new_maxA - new_minA) + new_minA
maxEnrgImg = max(max(enrgImg));
minEnrgImg = min(min(enrgImg));
enrgImg = (enrgImg - minEnrgImg)/(maxEnrgImg - minEnrgImg) * (0 - (-1)) + (-1);

% Show the image energy (FOR TESTING PURPOSE)
%figure, imshow(enrgImg,[]), title('Image energy')
%pause(1)
% Find the derivatives of the image energy function
fy = fspecial('sobel');
fx = fy';
gradY = imfilter(enrgImg, fy, 'replicate');
gradX = imfilter(enrgImg, fx, 'replicate');

% Initialize the alpha and beta values for each snake point
points(3,:) = alpha;
points(4,:) = beta;
h = 1;

% Construct the A matrix
A = constructA(points, n, h);
invAI = inv(A + (1/delta)*eye(n));

% Evolve the snake
for i = 1:maxIt
    if strcmpi(rs, 'on')
        % Resample snake for each 15 iterations
        if mod(i, 15) == 0
            [points] = snakeResample([[x y]'; points(3:4,:)]);
            n = size(points, 2);
            x = points(1,:)';
            y = points(2,:)';
            A = constructA(points, n, h);
            invAI = inv(A + (1/delta)*eye(n));
        end
        scale = 0.47;
    end
    % Interpolate to find the image energy on the discrete grid
    enrgX = interp2(gradX, x, y, '*linear');
    enrgY = interp2(gradY, x, y, '*linear');
    % Calculate new snake point indices
    x = invAI * ((1/delta)*x + enrgX);
    y = invAI * ((1/delta)*y + enrgY);
    % Stopping criterion check
    if i > 1
        previous = cellArray{i-1};
        if size(previous(1:2,:)) == size([x y]')
            if (norm([x y]' - previous(1:2,:))/n)*scale < stopThresh
                break;
            end
        end
    end
    cellArray(i) = {[[x y]'; points(3:4,:)]};
    stop = i;
end
fprintf('Number of snake control points at termination: %d \n', n);
% Delete cells that were not used
cellArray(stop+1:end) = [];
end

%% constructA
%
%
% Arguments:
%
% % % % % % %
% Subfunction that assembles the coefficient matrix A.
% points - A matrix of vectors, where the first row of each
% column contains the x coordinate of the point and the
% second row of each column contains the y coordinate of the
% point.
% n - The total number of control points along the snake
% curve.
% h - Small constant used in the finite difference
% approximation.

% Output:        A - The coefficient matrix A.
function [A] = constructA(points, n, h)
    
alphaP1 = [points(3,n) points(3,1:n-1)];
betaM1 = [points(4,2:n) points(4,1)];
betaP1 = [points(4,n) points(4,1:n-1)];

% Produce the five diagonal vectors
a = betaM1/h^4;
b = -2*(points(4,:)+betaM1)/h^4 - points(3,:)/h^2;
c = (betaP1+4*points(4,:)+betaM1)/h^4 + (alphaP1+points(3,:))/h^2;
d = -2*(betaP1+points(4,:))/h^4 - alphaP1/h^2;
e = betaP1/h^4;

% Assemble the A matrix
B = [d' e' a' b' c' d' e' a' b'];
B = fliplr(B);
A = full(spdiags(B, [-(n-1) -(n-2) -2 -1 0 1 2 (n-2) (n-1)], n, n))';
end


