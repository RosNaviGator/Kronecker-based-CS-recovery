function [X, nbVect, LstCPT, longs] = generateWaveletDictionary(DIM, varargin)
    %% Parameters
    LstCPT_DEF = {{'sym4', 5}};  % Default to Symlet 4 with 5 levels if not specified

    % Allowed wavelet names
    allowedWavelets = {'haar', 'db1', 'db2', 'db3', 'db4', 'db5', 'db6', ...
                       'db7', 'db8', 'db9', 'db10', 'sym2', 'sym3', 'sym4', ...
                       'sym5', 'sym6', 'sym7', 'sym8', 'coif1', 'coif2', ...
                       'coif3', 'coif4', 'coif5', 'bior1.1', 'bior1.3', ...
                       'bior1.5', 'bior2.2', 'bior2.4', 'bior2.6', 'bior2.8', ...
                       'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7', 'bior3.9', ...
                       'bior4.4', 'bior5.5', 'bior6.8', 'rbio1.1', 'rbio1.3', ...
                       'rbio1.5', 'rbio2.2', 'rbio2.4', 'rbio2.6', 'rbio2.8', ...
                       'rbio3.1', 'rbio3.3', 'rbio3.5', 'rbio3.7', 'rbio3.9', ...
                       'rbio4.4', 'rbio5.5', 'rbio6.8', 'dmey'};
    
    % Allowed extension modes
    allowedExtModes = {'zpd', 'sym', 'symw', 'asym', 'asymw', 'sp0', 'sp1', ...
                       'ppd', 'per'};
    
    % Handle input arguments
    if nargin > 1
        % Convert any string input to char arrays if needed
        [varargin{:}] = wavelet.internal.wconvertStringsToChars(varargin{:});
    end

    % Parse input arguments
    nbIN = length(varargin);
    k = 1;
    while k <= nbIN
        argNAM = lower(varargin{k});
        switch argNAM
            case 'lstcpt'
                LstCPT = varargin{k + 1}; 
                k = k + 2;
            otherwise
                error('Unknown argument name.');
        end
    end

    % Use defaults if not provided
    if ~exist('LstCPT', 'var'), LstCPT = LstCPT_DEF; end

    % Validate the wavelet settings
    for i = 1:length(LstCPT)
        cptname = LstCPT{i};
        if iscell(cptname)
            waveletName = cptname{1};
            if ~ismember(waveletName, allowedWavelets)
                error(['Invalid wavelet name: ', waveletName]);
            end
            
            if length(cptname) > 2
                extMode = cptname{3};
                if ~ismember(extMode, allowedExtModes)
                    error(['Invalid extension mode: ', extMode]);
                end
            end
        else
            error("Each entry in LstCPT must be a cell array containing at least the wavelet name and the level. E.g. {{'haar', 10}}");
        end
    end

    %% Dictionary Generation
    nbCOMPO = length(LstCPT);
    X = [];
    nbVect = [];
    longs = cell(1, nbCOMPO);

    for kCOMPO = 1:nbCOMPO
        cptname = LstCPT{kCOMPO};

        % Handle wavelet basis
        if iscell(cptname)
            level = cptname{2};
            if length(cptname) > 2
                dwtEXTM = cptname{3};
            end
            cptname = cptname{1}; 
        end

        % Set DWT extension mode if necessary
        if ~exist('dwtEXTM', 'var')
            dwtEXTM = 'per';
        end
        old_dwtEXTM = dwtmode('status', 'nodisp');
        dwtmode(dwtEXTM, 'nodisp');

        % Set level of decomposition if not specified
        if ~exist('level', 'var')
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

        % Restore the DWT extension mode if necessary
        dwtmode(old_dwtEXTM, 'nodisp');

        % Normalize the basis vectors
        S = sum(Z .* Z, 1);
        Z = Z ./ repmat(S.^0.5, DIM, 1);
        X = [X Z]; %#ok<AGROW>
        nbVect = [nbVect size(Z, 2)]; %#ok<AGROW>
    end

    % Convert X to sparse format
    X = sparse(X);
end
