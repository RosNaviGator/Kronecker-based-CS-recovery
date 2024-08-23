%This code shows the effect of Kronecker technique compressive sensing
%recovery; this technique has been used for different signals and measurement matrices, 
% and it has been published in multiple journals and conference papers:

% [1] H. Zanddizari, S. Rajan, and H. Zarrabi, “Increasing the quality of reconstructed  signal  in  compressive  sensing  utilizing  Kronecker  technique,”
     % Biomedical Engineering Letters, vol. 8, no. 2, pp. 239–247, May 2018.

% [2] D. Mitra, H. Zanddizari, and S. Rajan, "Investigation of Kronecker-based recovery of compressed ECG signal," 
     % IEEE Transactions on Instrumentation and Measurement, pp. 1-1, 2019.

% [3] D. Mitra, H. Zanddizari, and S. Rajan, “Improvement of signal quality during recovery of compressively sensed ECG signals,” 
     % in 2018 IEEE International Symposium on Medical Measurements and Applications (MeMeA), June 2018, pp. 1–5.

% [4] D. Mitra, H. Zanddizari, and S. Rajan, “Improvement of recovery in segmentation-based parallel compressive sensing,” 
     % in 2018 IEEE International Symposium on Signal Processing and Information Technology (ISSPIT), Dec 2018, pp. 501–506.

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


clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%
data0 = load('100m.mat');%ECG signal from MIT database
signal = data0.val(1,512*0+1:512*4)';%some samples of an ECG signal
%------------------------measurement matrix
CR = 1/4;%compression ratio 
N1 = 16; %length of signal for sensing 
M1= N1 * CR; %length of measurement 

N2 = 512; %length of signal for recovery based on Kronecker technique 
M2 = N2*CR;%length of concatened measurements
KronFact = N2/N1; % kronecker factor 

%------------------------------sensing matrix
%----random matrices
% A=((1/M2)^2)*randn(M2,N2);
%A=ones(M2,N2);A=binornd(A,.5);A=A-.5;A=1/sqrt(M2)*A;
%------ Deterministic matrix: DBBD matrix 
A = zeros(M1,N1);m=N1/M1;for i=1:M1 A(i,1+(i-1)*m:(i)*m) = 1;end
%------------------------------Sparsifying basis matrix
dic = 1;
if dic == 1
dict1 = wmpdictionary(N1,'LstCpt',{'dct'});
dict2 = wmpdictionary(N2,'LstCpt',{'dct'});

elseif dic==0
dict1 = wmpdictionary(N1,'lstcpt',{{'db10',10}});
dict2 = wmpdictionary(N2,'lstcpt',{{'db10',10}});
end

%----------Generating large sesing matrix for recovery, not for sensing
AA = kron(eye(KronFact),A);
A1 = A*dict1;
A2 = AA*dict2;

%------------------------Compression phase based of CS
for i = 1:length(signal)/N1
   Y(:,i)=A*signal((i-1)*N1+1:N1*i,1);
end

%------------------------initializing Sl0 parameter
for i=1:1
    sigma_off = 0.001;
    A_pinv1 = pinv(A1);A_pinv2 = pinv(A2);mu_0 = 2;sigma_decrease_factor = 0.5;L = 3;
    true_s = sparseSigGen4plusNoise(9,floor(27/4),sigma_off);
    if sigma_off>0
        sigma_min = sigma_off*4;
    else
        sigma_min = 0.00001;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%Recovery process%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------- Ordinary CS recovery method 
for i=1:length(signal)/N1
        y = Y(:,i);
        xp = SL0(A1, y, sigma_min, sigma_decrease_factor, mu_0, L, A_pinv1);
        %xp = l1eq_pd(A'*y2, A2, [], y2, 1e-1);
        Recovered_ordinary(N1*i-(N1-1):N1*i) = dict1*xp;
end
%------------------------ Kroneckered-based recovery
for i=1:length(signal)/N2
        y_concatened = reshape(Y(:,(i-1)*KronFact+1:i*KronFact),[M2,1]);
        % y_concatened: contain multiple measurement vectors without any change
        xp=SL0(A2, y_concatened, sigma_min, sigma_decrease_factor, mu_0, L, A_pinv2);
        %xp = l1eq_pd(AA'*y1, A1, [], y1, 1e-1);
        Recovered_Kronecker(N2*i-(N2-1):N2*i) = dict2*xp;
end

% %-----------------------------------------Calculating SNR
% disp('SNR of ordinary CS recovery');
err1 = Recovered_ordinary-signal';%new sensing
SNR_Ordinary= 20*log10(norm(signal)/norm(err1))

% disp('SNR of Kronecker-based recovery');
err2 = Recovered_Kronecker-signal';
SNR_Proposed= 20*log10(norm(signal)/norm(err2))


figure;
subplot(2,1,1);plot(signal);hold on;plot(Recovered_ordinary,'r');title('Ordinary CS recovery');
subplot(2,1,2);plot(signal);hold on;plot(Recovered_Kronecker,'k');title('Kronecker-based recovery');
