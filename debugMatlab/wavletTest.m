clc, clear, close all


%% PARAMETERS
% Initialize the parameters based on the function call
DIM = 4;  % Example dimension, this should match your signal length or desired basis size
LstCPT = {{'haar', 10}};  % The list of components, one used in the project LstCPT = {{'db10', 10}};
printThem = true;  % Print the matrices if they are small enough



%% USE LIBRARY DIRECTLY
% Create a dictionary
dwt = wmpdictionary(DIM, 'lstcpt', LstCPT);
% Make non-sparse
dwt = full(dwt);


%% LIBRARY IMITATION FUNCTION
% Create a dictionary
my_dwt = generateWaveletDictionary(DIM,'lstcpt', LstCPT);
% Make non-sparse
my_dwt = full(my_dwt);


%% LIBRARY IMITATION SCRIPT
% Variables to start the code snippet you provided
nbCOMPO = length(LstCPT);
X = [];
nbVect = [];
longs = cell(1, nbCOMPO);


for kCOMPO = 1:nbCOMPO
    cptname = LstCPT{kCOMPO};
    
    % Wavelet basis (excluding wavelet packet)
    if iscell(cptname)
        level = cptname{2};
        if length(cptname) > 2
            dwtEXTM = cptname{3};
        end
        cptname = cptname{1};
    end
    
    % Change the DWT extension mode if necessary.
    if ~exist('dwtEXTM','var')
        dwtEXTM = 'per';
    end
    old_dwtEXTM = dwtmode('status', 'nodisp');
    dwtmode(dwtEXTM, 'nodisp');
    
    % Get the level of decomposition.
    if ~exist('level','var')
        level = 5;
    end
    
    % Perform wavelet decomposition
    [~, lon] = wavedec(zeros(DIM, 1), level, cptname);
    NbCFS = sum(lon(1:end-1));
    DEC = mdwtdec('c', zeros(DIM, NbCFS), level, cptname);
    sCa = size(DEC.ca);
    for k = 1:sCa(1)
        DEC.ca(k, k) = 1;
    end
    icol = sCa(1) + 1;
    for j = level:-1:1
        sCd = size(DEC.cd{j});
        for k = 1:sCd(1)
            DEC.cd{j}(k, icol) = 1;
            icol = icol + 1;
        end
    end
    longs{kCOMPO} = lon;
    Z = mdwtrec(DEC);
    
    clear dwtEXTM level
    
    % Restore the DWT extension mode if necessary.
    dwtmode(old_dwtEXTM, 'nodisp');
    
    % Normalize the basis vectors (applies to both wavelet and others)
    S = sum(Z .* Z, 1);
    Z = Z ./ repmat(S.^0.5, DIM, 1);
    X = [X Z];               %#ok<AGROW>
    nbVect = [nbVect size(Z, 2)];   %#ok<AGROW>
end

% X = sparse(X);






%% EVALUATION
% control if function is equal to original entry by entry
if (size(dwt) == size(my_dwt))
    if (all(dwt(:) == my_dwt(:)))
        disp('Function and original are equal')
    else
        disp('Function and original are not equal')
    end
else
    disp('Function and original are not equal')
end

% control if script is equal to original entry by entry
if (size(dwt) == size(X))
    if (all(dwt(:) == X(:)))
        disp('Script and original are equal')
    else
        disp('Script and original are not equal')
    end
else
    disp('Script and original are not equal')
end

% control if the two imitated are equal
if (size(my_dwt) == size(X))
    if (all(my_dwt(:) == X(:)))
        disp('Function and script are equal')
    else
        disp('Function and script are not equal')
    end
else
    disp('Function and script are not equal')
end

fprintf('\n')


fprintf('From library (original wmptdictionary):\n')
disp(['DIM: ', num2str(DIM), '| LstCPT: ', char(LstCPT{1}{1}), '| Level: ', num2str(LstCPT{1}{2}), '| Size of dwt: ', num2str(size(dwt))]);
if (DIM <= 32 && printThem)
    % Print dwt entirely
    disp(dwt)
end

fprintf('Imitated with function:\n')
disp(['DIM: ', num2str(DIM), '| LstCPT: ', char(LstCPT{1}{1}), '| Level: ', num2str(LstCPT{1}{2}), '| Size of my_dwt: ', num2str(size(my_dwt))]);
if (DIM <= 32 && printThem)
    % Print my_dwt entirely
    disp(my_dwt)
end

fprintf('Imitated with script:\n')
disp(['DIM: ', num2str(DIM), '| LstCPT: ', char(LstCPT{1}{1}), '| Level: ', num2str(LstCPT{1}{2}), '| Size of X: ', num2str(size(X))]);
if (DIM <= 32 && printThem)
    % print X entirely
    disp(X)
end

