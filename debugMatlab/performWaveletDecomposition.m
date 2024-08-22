%% Unravel dictionary creation

clc, clear, close all


%% Perform wavelet decomposition
DIM = 4;
cptname = 'db6';
level = 7;
dwtmode('per'); 


[~, lon] = wavedec(zeros(DIM, 1), level, cptname);
disp(lon)
% lon is a vector of the number of coefficients at each level of the decomposition


NbCFS = sum(lon(1:end-1));
disp(NbCFS)



DEC = mdwtdec('c', ones(DIM, NbCFS), level, cptname);
disp(DEC)
% 1. Print the decomposition direction
disp('Decomposition Direction:');
disp(DEC.dirDec);
% 2. Print the decomposition level
disp('Decomposition Level:');
disp(DEC.level);
% 3. Print the wavelet name
disp('Wavelet Name:');
disp(DEC.wname);
% 4. Print the dwtFilters structure
disp('DWT Filters:');
disp(DEC.dwtFilters);
% 5. Print the signal extension mode
disp('DWT Extension Mode:');
disp(DEC.dwtEXTM);
% 6. Print the dwtShift
disp('DWT Shift:');
disp(DEC.dwtShift);
% 7. Print the size of the data
disp('Data Size:');
disp(DEC.dataSize);
% 8. Print the approximation coefficients (ca)
disp('Approximation Coefficients (ca):');
disp(DEC.ca);
% 9. Print the detail coefficients (cd)
disp('Detail Coefficients (cd):');
for i = 1:length(DEC.cd)
    fprintf('Detail Coefficients at Level %d:\n', i);
    disp(DEC.cd{i});
end



sCa = size(DEC.ca);
fprintf('sCa: ');
disp(sCa)


fprintf('DEC.ca: ');
disp(DEC.ca)
for k = 1:sCa(1)
    DEC.ca(k, k) = 1;
end
fprintf('DEC.ca: ');
disp(DEC.ca)






icol = sCa(1) + 1;
for j = level:-1:1
    sCd = size(DEC.cd{j});
    for k = 1:sCd(1)
        DEC.cd{j}(k, icol) = 1;
        icol = icol + 1;
    end
end

%longs{kCOMPO} = lon;


Z = mdwtrec(DEC);
disp('Z:');
disp(Z)


% Normalize the basis vectors
S = sum(Z .* Z, 1);
Z = Z ./ repmat(S.^0.5, DIM, 1);
disp('Normalized Z:');
disp(Z)

