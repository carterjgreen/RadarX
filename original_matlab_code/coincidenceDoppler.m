function [PosSpeeds, NegSpeeds] = coincidenceDoppler(BinCell,Nmax)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Interpret Doppler Bins
n_PRF = size(BinCell, 2);
n_targets = sum(nnz(BinCell{1}));
%% Find Bin Lengths and LCM
L = binLength(BinCell);
Lu_max = lcms(L); % Downloaded from Mathworks library
%% How many targets in each row
n_det = zeros(1, size(BinCell, 2));
for i = 1:size(BinCell, 2)
    n_det(i) = nnz(BinCell{i});
end

%% Detect Zero Points
ZeroPoint = cellfun(@length, BinCell)/2 + 1;

%% Make new Positive Structure
for i = 1:size(BinCell,2)
    n_rep = floor(Nmax/length(BinCell{i})) + 1; % Quotient + 1
    PosCell{i} = repmat(BinCell{i}, [1 n_rep]);
end
%% Take Out First Half
for i = 1:length(PosCell)
    PosCell{i}(1:ZeroPoint(i)-1) = [];
end
%% Check Detections Positive
zeroFlag = zeroTest(BinCell, n_PRF);
% Flag 1
flag1 = ~all(n_det == n_det(1));
if ~zeroFlag 
    detections = detection(PosCell, n_PRF);
    PosSpeeds = binMatch2(n_PRF, detections);
    % Flag 2
%     aTargets = PosSpeeds(Lu_max < PosSpeeds & PosSpeeds< Nmax);
%     flag2 = any(aTargets);
    % Flag 3
%     flag3 = length(PosSpeeds) > max(n_det);
else
    PosSpeeds = [];
end

%% Create Negative Bins
for i = 1:length(BinCell)
    BinCell{i} = fliplr(BinCell{i}); 
end
%% Make new Negative Structure
for i = 1:size(BinCell,2)
    n_rep = floor(Nmax/length(BinCell{i})) + 1; % Quotient + 1
    NegCell{i} = repmat(BinCell{i}, [1 n_rep]);   
end
%% Remove First Half
for i = 1:length(NegCell)
   NegCell{i}(1:ZeroPoint(i)) = []; 
end
%% Check Detections Negative
zeroFlag = zeroTest(BinCell, n_PRF);
% Flag 1
flag1 = ~all(n_det == n_det(1));
if ~zeroFlag 
    detections = detection(NegCell, n_PRF);
    NegSpeeds = binMatch2(n_PRF, detections);
    % Flag 2
%     aTargets = PosSpeeds(Lu_max < PosSpeeds & PosSpeeds< Nmax);
%     flag2 = any(aTargets);
    % Flag 3
%     flag3 = length(PosSpeeds) > max(n_det);
else
    NegSpeeds = [];
end

end

% Local Functions Start Here
function [L] = binLength(B)
%UNTITLED2 Summary of this function goes here
%   Finds length of each bin
n_PRF = size(B, 2);
%% Find Bin Lengths and LCM
L = zeros(1, n_PRF);
for i = 1:n_PRF
    L(i) = size(B{i}, 2);
end
end

function [flag] = zeroTest(B, n_PRF)
%   Checks to see if a PRF has no returns while others do
flag_array = zeros(1, n_PRF);
for i = 1:n_PRF
    if any(B{i})
        flag_array(i) = 0;
    else
        flag_array(i) = 1;
    end
end
flag = any(flag_array);
end

function output = lcms(numberArray)

numberArray = reshape(numberArray, 1, []);

%% prime factorization array
for i = 1:size(numberArray,2)
    temp = factor(numberArray(i));
   
    for j = 1:size(temp,2)
        output(i,j) = temp(1,j);
    end
end

%% generate prime number list
p = primes(max(max(output)));
%% prepare list of occurences of each prime number
q = zeros(size(p));

%% generate the list of the maximum occurences of each prime number
for i = 1:size(p,2)
    for j = 1:size(output,1)
        temp = length(find(output(j,:) == p(i)));
        if(temp > q(1,i))
            q(1,i) = temp;
        end
    end
end

%% the algorithm
z = p.^q;

output = 1;

for i = 1:size(z,2)
    output = output*z(1,i);
end
end

function [cell] = detection(B, n_PRF)
%UNTITLED3 Summary of this function goes here
%   Returns a cell that contains a list of detections in each bin 
cell = {};
for i = 1:n_PRF
  cell{i} = find(B{i} == 1);  
end
end

function [match] = binMatch2(n_PRF, detectionCell)
%UNTITLED6 Summary of this function goes here
%   This function finds the locations of intersections in the PRFs
match = intersect(detectionCell{1}, detectionCell{2});
if(size(detectionCell, 2) > 2)
    for i = 3:n_PRF
        match = intersect(match, detectionCell{i});
    end
end

end


