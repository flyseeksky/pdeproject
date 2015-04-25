
function main(type, input, ssc)
close all

% Set image path
pathName = './';

% Define name of image file
imgName = 'bluestar.jpg';

% Load the image file
originalImg = imread([pathName imgName]);

% Convert image to grayScale if necessary
fileInfo = imfinfo([pathName imgName]);
if strcmp(fileInfo.ColorType, 'truecolor')
    img = rgb2gray(originalImg);
else
    img = originalImg;
end

% If input = user, then let the user input the snake control points
if strcmpi(input, 'user')
    
    % Create the input figure window
    figure('Name', 'Snake: Insert control points', 'NumberTitle', 'off', 'Resize', 'on');
    imshow(originalImg);
    title('\fontsize{13}\fontname{Monospaced} Click outside of image to finish.'); hold on;
    % Counter for the total number of points
    numPoints = 0;
    % Input loop - lets user enter points
    while (true)
        % Wait for mouse input
        waitforbuttonpress;

        % Get the indices of the new point point
        newP=get(gca,'CurrentPoint');
        newP=newP(1,1:2);

        % Check if point is outside of the image
        if (newP(2)>size(img, 1) || newP(1)>size(img, 2) || ...
            newP(2)< 1 || newP(1)< 1)
            break; 
        end
        numPoints = numPoints + 1;
        points(1, numPoints)=(newP(1));
        points(2, numPoints)=(newP(2));

        % Plot the new control point
        plot(points(1, numPoints),points(2, numPoints), ...
            'o', ...
            'MarkerSize', 5, ...
            'MarkerFaceColor', 'b', ...
            'MarkerEdgeColor', 'b');
   end

    %points = load('points.txt');
    %savefile = 'points.txt';
    %save(savefile, 'points', '-ASCII');
else
    % Plot the snake as a circle around the object in the image
    % Circle radius
    r = 70;
    % Cicle center [x y]
    c = [150 110];
    % Number of points in the snake
    N = 50;
    % Calculate snake points in a circle
    points(1, :) = c(1)+floor(r*cos((1:N)*2*pi/N)+0.5);
    points(2, :) = c(2)+floor(r*sin((1:N)*2*pi/N)+0.5);
end

C = {points};

% Show initial snake
showSnake(C, originalImg);
pause(0.1);
figure;

% Run the snake algorithm
if strcmpi(type, 'Kass')
    % alpha: controls elasticity
    % beta:  controls curvature
    % delta: controls step size
    % sigma: controls amount of Gaussian blurring
    % maxIt: Defines the maximum number of snake iterations
    % rs:    Controls whether to have resmapling on or off
    
    alpha = 0.05; beta = 0.0005; delta = 1; sigma = 3; maxIt = 16000; rs = 'on';
    
    %alpha = 0.035;
    % Use scale space continuation if ssc = on
    if strcmpi(ssc, 'on')
        sigma = 15;
        iteration = ceil(maxIt/4);
        C = cell(0);
        for i = 1:4
            Ctmp = KassSnake(points, img, alpha, beta, delta, sigma, iteration, rs);
            points = Ctmp{end};
            sigma = sigma - 4;
            C = [C; Ctmp];
        end
    else
        C = KassSnake(points, img, alpha, beta, delta, sigma, maxIt, rs);
    end
    fprintf('\n');
    fprintf('Running %s snake with following parameters: \n', type);
    fprintf('Alpha: %d \n', alpha);
    fprintf('Beta: %d \n', beta);
    fprintf('Delta: %d \n', delta);
    fprintf('Sigma: %d \n', sigma);
    fprintf('Max iterations: %d \n', maxIt);
    fprintf('Resampling: %s \n', rs);
    fprintf('Scale space continuation: %s \n', ssc);
elseif strcmpi(type, 'greedy')
    % alpha: controls continuity (higher numbers result in points being
    %                             equally spaced.)
    % beta:  controls curvature
    % gamma: controls strength of image energy
    % s:     controls the size of the neighborhood
    % sigma: controls amount of Gaussian blurring
    % maxIt: Defines the maximum number of snake iterations
    alpha = 1.2; beta = 1; gamma = 1.2; s = 5; sigma = 3; maxIt = 200;
    % Use scale space continuation if ssc = on
    if strcmpi(ssc, 'on')
        sigma = 15;
        iteration = ceil(maxIt/4);
        C = cell(0);
        for i = 1:4
            Ctmp = GreedySnake(points, img, alpha, beta, gamma, s, sigma, iteration);
            points = Ctmp{end};
            sigma = sigma - 4;
            C = [C; Ctmp];
        end
    else
        C = GreedySnake(points, img, alpha, beta, gamma, s, sigma, maxIt);
    end
    fprintf('\n');
    fprintf('Running %s snake with following parameters: \n', type);
    fprintf('Alpha: %d \n', alpha);
    fprintf('Beta: %d \n', beta);
    fprintf('Gamma: %d \n', gamma);
    fprintf('Neighborhood size: %d \n', s);
    fprintf('Sigma: %d \n', sigma);
    fprintf('Max iterations: %d \n', maxIt);
    fprintf('Scale space continuation: %s \n', ssc);
end
% Show snake evolution
figure
showSnake(C, originalImg, beta);
clear all
end




function showSnake(C, img, betaStart)
it = size(C, 1);
if it < 200
    step = 1;
else
    step = floor(it/100);
end
for i = 0:step:it
    cla;
    imshow(img);
    hold on;
    if it == 1
        title('\fontsize{13}\fontname{Monospaced} Initial snake');
        points = C{1};
    elseif i == 0
        points = C{i+1};
        title(['\fontsize{13}\fontname{Monospaced} Iteration number: ', num2str(i+1)]);
    else
        points = C{i};
        title(['\fontsize{13}\fontname{Monospaced} Iteration number: ', num2str(i)]); 
    end
end
    
% Plot lines between points
plot(points(1,:), points(2,:), ...
    '-g', ...
    'LineWidth', 3);
plot([points(1,end) points(1,1)], [points(2,end) points(2,1)], ...
    '-g', ...
    'LineWidth', 3);

% Plot points, change colour of points where
% the beta value is relaxed
for j = 1:size(points, 2);
    if size(points, 1) == 2
        plot(points(1,j), points(2,j), ...
            'o', ...
            'MarkerSize', 7, ...
                            'MarkerFaceColor', 'b', ...
            'MarkerEdgeColor', 'b');
    elseif points(4,j) ~= betaStart
        plot(points(1,j), points(2,j), ...
            'o', ...
            'MarkerSize', 7, ...
            'MarkerFaceColor', 'r', ...
            'MarkerEdgeColor', 'r');
    else
        plot(points(1,j), points(2,j), ...
            'o', ...
            'MarkerSize', 7, ...
            'MarkerFaceColor', 'b', ...
            'MarkerEdgeColor', 'b');
    end
end

% Delay between each iteration
pause(0.005);
end