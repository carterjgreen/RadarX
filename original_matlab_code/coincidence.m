function [Targets] = coincidence(BinCell,Nmax)
%UNTITLED Applies Coincidence Algorithm to Diasmbiguate targets in range or
% velocity
%Ts = 1e-6;
%% Interpret Range Bins
n_PRF = size(BinCell, 2);
n_targets = sum(nnz(BinCell{1}));
%% Find Bin Lengths and LCM
L = binLength(BinCell);
Lu_max = lcms(L); % Downloaded from Mathworks library
% PRF1 = 1/(length(B{1})*Ts);
%% How many targets in each row
n_det = zeros(1, size(BinCell, 2));
for i = 1:size(BinCell, 2)
    n_det(i) = nnz(BinCell{i});
end
%% Make new structure
for i = 1:size(BinCell,2)
    n_rep = floor(Nmax/length(BinCell{i})) + 1; % Quotient + 1
    myCell{i} = repmat(BinCell{i}, [1 n_rep]);
end
%% Check Detections
zeroFlag = zeroTest(BinCell, n_PRF);
%% Flag 1
flag1 = ~all(n_det == n_det(1));
if ~zeroFlag 
    detections = detection(myCell, n_PRF);
    Targets = binMatch2(n_PRF, detections);
%     targetDist = delR*Targets;
    %% Flag 2
    aTargets = Targets(Lu_max < Targets & Targets< Nmax);
    flag2 = any(aTargets);
    %% Flag 3
    flag3 = length(Targets) > max(n_det);
    %% Flag 4
    % flag4 = length
else
    Targets = [];
    %flags = []
end

%flags = [flag1 flag2 flag3];
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

