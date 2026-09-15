% Coordinates of multiple objects of interest
objectsCoords = [];

% Coordinates of neighbors (same for all objects)
neighborsCoords = [];
%%

% Calculate distances for all objects to all neighbors
distances = calculateDistances(objectsCoords, neighborsCoords);

% Find the minimum distance for each object
[minDistances, idx] = min(distances, [], 2);

% Prepare results array
results = zeros(size(objectsCoords, 1), 5);  % Initialize results array with correct number of columns

% Store results
for i = 1:size(objectsCoords, 1)
    nearestNeighborCoords = neighborsCoords(idx(i), :);
    results(i, :) = [objectsCoords(i, 1), objectsCoords(i, 2), nearestNeighborCoords(1), nearestNeighborCoords(2), minDistances(i)];
end

% Display results array
disp(results);

% Function to calculate Euclidean distances for multiple objects
function distances = calculateDistances(objectsCoords, neighborsCoords)
    numObjects = size(objectsCoords, 1);
    numNeighbors = size(neighborsCoords, 1);
    
    % Initialize the distance matrix
    distances = zeros(numObjects, numNeighbors);
    
    % Calculate distances for each object to all neighbors
    for i = 1:numObjects
        % Differences for current object
        differences = neighborsCoords - objectsCoords(i, :);
        
        % Squared distances
        squaredDistances = differences .^ 2;
        
        % Sum and square root to get Euclidean distances
        distances(i, :) = sqrt(sum(squaredDistances, 2));
    end
end