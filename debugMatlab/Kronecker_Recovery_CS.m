%This code shows the effect of Kronecker technique compressive sensing
%recovery; this technique has been used for different signals and measurement matrices, 
% and it has been published in multiple journals and conference papers:

% [1] H. Zanddizari, S. Rajan, and H. Zarrabi, �Increasing the quality of reconstructed  signal  in  compressive  sensing  utilizing  Kronecker  technique,�
     % Biomedical Engineering Letters, vol. 8, no. 2, pp. 239�247, May 2018.

% [2] D. Mitra, H. Zanddizari, and S. Rajan, "Investigation of Kronecker-based recovery of compressed ECG signal," 
     % IEEE Transactions on Instrumentation and Measurement, pp. 1-1, 2019.

% [3] D. Mitra, H. Zanddizari, and S. Rajan, �Improvement of signal quality during recovery of compressively sensed ECG signals,� 
     % in 2018 IEEE International Symposium on Medical Measurements and Applications (MeMeA), June 2018, pp. 1�5.

% [4] D. Mitra, H. Zanddizari, and S. Rajan, �Improvement of recovery in segmentation-based parallel compressive sensing,� 
     % in 2018 IEEE International Symposium on Signal Processing and Information Technology (ISSPIT), Dec 2018, pp. 501�506.

%Author Hadi Zanddizari, 
% hadiz@mail.usf.edu
% hadizand@alumni.iust.ac.ir


%-----------The main objective of this approach
% For fast and efficient compression, sensing phase in compressive sensing
%can be done in very small size,
% because in this case: 
    %it requires very small measurement matrix, 
    %less number of multiplication and addition operations, and
    %less delay for generating compressed samples
    % which can be used for sensors with low computational resources.
%But sensing in small size degrades quality of recovery. Kronecker
%technique can be used in order to improve the quality of recovered signal.

%---------------------------------------CS Recovery algorithm
%Any recovery algorithm can be used. In this code, Sl0 which is very fast CS recovery algorithm is used.
%[sl0-reference]: http://ee.sharif.edu/~SLzero/

%----------------------------------Database
%In this code, an ECG signal from MIT-BIH Arrhythmia database which is public dataset is used.
%[Database Reference]:MIT-BIH Arrhythmia Database. [Online]. Available: http://www.physionet.org/physiobank/database/mitdb/
%--------------------------------------------------------


clear ; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%
data0 = load('100m.mat');%ECG signal from MIT database
signal = data0.val(1,512*0+1:512*4)';%some samples of an ECG signal

% print signal shape
disp(['Size of signal: ', num2str(size(signal))]);

%------------------------measurement matrix
CR = 1/4;%compression ratio 
N1 = 16; %length of signal for sensing 
M1= N1 * CR; %length of measurement 

N2 = 512; %length of signal for recovery based on Kronecker technique 
M2 = N2*CR;%length of concatened measurements
KronFact = N2/N1; % kronecker factor 

%------------------------------sensing matrix
%----random matrices

RAND = 0;
if RAND == 1
    A = ((1/M1)^2)*randn(M1,N1);
end


BINOMIAL = 0;
if BINOMIAL == 1    
    A = ones(M1,N1);
    A=binornd(A,.5);
    A=A-.5;
    A=1/sqrt(M1)*A;
end

%------ Deterministic matrix: DBBD matrix 
DBBD = 1;
if DBBD == 1
A = zeros(M1,N1);
    m=N1/M1;
    for i=1:M1 
        A(i,1+(i-1)*m:(i)*m) = 1;
    end

end

% A random selection of -1, +1
RANDONES = 0;
if RANDONES == 1
    A = 2*randi([0 1],M1,N1)-1;
end



%------------------------------Sparsifying basis matrix
dic = 1;
if dic == 1
dict1 = wmpdictionary(N1,'LstCpt',{'dct'});
dict2 = wmpdictionary(N2,'LstCpt',{'dct'});

elseif dic==0
dict1 = wmpdictionary(N1,'lstcpt',{{'db10',10}});
dict2 = wmpdictionary(N2,'lstcpt',{{'db10',10}});
end





%----------Generating large sensing matrix for recovery, not for sensing
AA = kron(eye(KronFact),A);
A1 = A*dict1;
A2 = AA*dict2;



%------------------------Compression phase based of CS
Y = zeros(M1,length(signal)/N1);
for i = 1:length(signal)/N1
   Y(:,i)=A*signal((i-1)*N1+1:N1*i,1);
end



%------------------------initializing Sl0 parameter
for i=1:1
    
    sigma_off = 0.001;
    A_pinv1 = pinv(A1);
    A_pinv2 = pinv(A2);
    mu_0 = 2;
    sigma_decrease_factor = 0.5;
    L = 3;
    true_s = sparseSigGen4plusNoise(9,floor(27/4),sigma_off);
    
    if sigma_off>0
        sigma_min = sigma_off*4;
    else
        sigma_min = 0.00001;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%Recovery process%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------- Ordinary CS recovery method 


% Create the directory if it doesn't exist
outputDir = 'debugCsvMAT';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% List of filenames to delete before starting
filesToDelete = {'y_block.csv', 's_block.csv', 'y_concatenated.csv', 's_block_kron.csv'};

% Delete the files if they exist
for i = 1:length(filesToDelete)
    filePath = fullfile(outputDir, filesToDelete{i});
    if exist(filePath, 'file')
        delete(filePath);
    end
end



% Define the precision for formatting
precision = '6';  % You can change '6' to any other number of decimal places


% Open the CSV files in append mode ('a' stands for append)
fileID_y = fopen(fullfile(outputDir, 'y_block.csv'), 'a');
fileID_xp = fopen(fullfile(outputDir, 's_block.csv'), 'a');

for i = 1:length(signal)/N1
    y = Y(:, i);
    xp = SL0(A1, y, sigma_min, sigma_decrease_factor, mu_0, L, A_pinv1);
    
    % Recover the signal
    Recovered_ordinary(N1*i-(N1-1):N1*i) = dict1*xp;
    
    % Generate the format specification string
    formatSpec_y = [repmat(['%', precision, 'f,'], 1, size(y, 1)-1), '%', precision, 'f\n'];
    formatSpec_xp = [repmat(['%', precision, 'f,'], 1, size(xp, 1)-1), '%', precision, 'f\n'];
    
    % Write y to the CSV file
    fprintf(fileID_y, formatSpec_y, y);
    
    % Write xp (s_block) to the CSV file
    fprintf(fileID_xp, formatSpec_xp, xp);
end

% Close the files
fclose(fileID_y);
fclose(fileID_xp);




%------------------------ Kroneckered-based recovery



% Define the precision for formatting
precision = '6';  % You can change '6' to any other number of decimal places

% Open the CSV files in append mode ('a' stands for append)
fileID_y_concat = fopen(fullfile(outputDir, 'y_concatenated.csv'), 'a');
fileID_xp_kron = fopen(fullfile(outputDir, 's_block_kron.csv'), 'a');

% Initialize the concatenated measurement vector
y_concatenated = zeros(M2, 1);

% Determine the length of the zero line
zeroLine_y = zeros(1, size(y_concatenated, 1));
zeroLine_xp = zeros(1, size(xp, 1));

for i = 1:length(signal)/N2
    % Reshape and concatenate the measurement vectors
    y_concatenated = reshape(Y(:, (i-1)*KronFact+1:i*KronFact), [M2, 1]);

    disp(['y_concatenated size: ', num2str(size(y_concatenated))]);
    disp(y_concatenated);
    
    % Perform SL0 algorithm (or any other recovery method)
    xp = SL0(A2, y_concatenated, sigma_min, sigma_decrease_factor, mu_0, L, A_pinv2);
    
    % Recover the signal
    Recovered_Kronecker(N2*i-(N2-1):N2*i) = dict2*xp;
    
    % Generate the format specification string for y_concatenated and xp
    formatSpec_y_concat = [repmat(['%', precision, 'f,'], 1, size(y_concatenated, 1)-1), '%', precision, 'f\n'];
    formatSpec_xp_kron = [repmat(['%', precision, 'f,'], 1, size(xp, 1)-1), '%', precision, 'f\n'];
    
    % Write y_concatenated to the CSV file
    fprintf(fileID_y_concat, formatSpec_y_concat, y_concatenated);
    
    % Write a zero line separator in y_concatenated.csv
    fprintf(fileID_y_concat, [repmat(['%', precision, 'f,'], 1, size(zeroLine_y, 2)-1), '%', precision, 'f\n'], zeroLine_y);
    
    % Write xp (s_block_kron) to the CSV file
    fprintf(fileID_xp_kron, formatSpec_xp_kron, xp);
    
    % Write a zero line separator in s_block_kron.csv
    fprintf(fileID_xp_kron, [repmat(['%', precision, 'f,'], 1, size(zeroLine_xp, 2)-1), '%', precision, 'f\n'], zeroLine_xp);
end

% Close the files
fclose(fileID_y_concat);
fclose(fileID_xp_kron);



% %-----------------------------------------Calculating SNR
% disp('SNR of ordinary CS recovery');
err1 = Recovered_ordinary-signal';%new sensing
SNR_Ordinary= 20*log10(norm(signal)/norm(err1));

% disp('SNR of Kronecker-based recovery');
err2 = Recovered_Kronecker-signal';
SNR_Proposed= 20*log10(norm(signal)/norm(err2));


figure;
subplot(2,1,1);plot(signal);hold on;plot(Recovered_ordinary,'r');title('Ordinary CS recovery');
subplot(2,1,2);plot(signal);hold on;plot(Recovered_Kronecker,'k');title('Kronecker-based recovery');

% print SNR
disp(['SNR of ordinary CS recovery: ',num2str(SNR_Ordinary)]);
disp(['SNR of Kronecker-based recovery: ',num2str(SNR_Proposed)]);


%------------------------End of the code------------------------%


% Create the directory if it doesn't exist
outputDir = 'debugCsvMAT';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Save signal to a CSV file with 6 decimal places
saveMatrixWithPrecision(signal, fullfile(outputDir, 'signal.csv'), '.6');

% Save A (Phi) to a CSV file with 6 decimal places
saveMatrixWithPrecision(A, fullfile(outputDir, 'A_Phi.csv'), '.6');

% Save AA (Phi_kron) to a CSV file with 6 decimal places
saveMatrixWithPrecision(AA, fullfile(outputDir, 'AA_Phi_kron.csv'), '.6');

% Save A1 (Theta) to a CSV file with 6 decimal places
saveMatrixWithPrecision(A1, fullfile(outputDir, 'A1_Theta.csv'), '.6');

% Save A2 (Theta_kron) to a CSV file with 6 decimal places
saveMatrixWithPrecision(A2, fullfile(outputDir, 'A2_Theta_kron.csv'), '.6');

% Save dict1 (Dict) to a CSV file with 6 decimal places
% Make non-sparse matrix
dict1 = full(dict1);
saveMatrixWithPrecision(dict1, fullfile(outputDir, 'dict1_Dict.csv'), '.6');

% Save dict2 (Dict_kron) to a CSV file with 6 decimal places
% Make non-sparse matrix
dict2 = full(dict2);
saveMatrixWithPrecision(dict2, fullfile(outputDir, 'dict2_Dict_kron.csv'), '.6');

% Save Y to a CSV file with 6 decimal places
saveMatrixWithPrecision(Y, fullfile(outputDir, 'Y.csv'), '.6');

% Save A_pinv1 to a CSV file with 6 decimal places
saveMatrixWithPrecision(A_pinv1, fullfile(outputDir, 'A_pinv1.csv'), '.6');

% Save A_pinv2 to a CSV file with 6 decimal places
saveMatrixWithPrecision(A_pinv2, fullfile(outputDir, 'A_pinv2.csv'), '.6');



% Function to save matrix with specific precision using fprintf
function saveMatrixWithPrecision(matrix, filename, precision)
    fileID = fopen(filename, 'w');
    formatSpec = [repmat(['%', precision, 'f,'], 1, size(matrix, 2)-1), '%', precision, 'f\n'];
    fprintf(fileID, formatSpec, matrix.');
    fclose(fileID);
end
