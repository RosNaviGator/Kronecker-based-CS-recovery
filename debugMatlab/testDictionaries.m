clear, clc, close all;

DIM = 4;

% Create a dictionary dct
dct = wmpdictionary(DIM,'LstCpt',{'dct'});
% make non sparse
dct = full(dct);

% Create a dictionary db10
dwt = wmpdictionary(DIM,'lstcpt',{{'db10',10}});
% make non sparse
dwt = full(dwt);

% print dct entirely
disp(['Size of dct: ', num2str(size(dct))]);
disp(dct)

% print dwt entirely
disp(['Size of dwt: ', num2str(size(dwt))]);
disp(dwt)

