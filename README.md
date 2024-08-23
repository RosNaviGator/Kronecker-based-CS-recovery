# Kronecker-based-CS-recovery -- ECG data -- RosNaviGator


## Fork to study and to port to python (if I manage) the work of [Kronecker-based-CS-recovery](https://github.com/hadizand/Kronecker-based-CS-recovery)


## Table of content of the README
- Original description (use it to understand what the repo is about, contains __references__ to papers about Kronecker Te)


## Original description (by [Hadi Zanddizari](https://github.com/hadizand))
Kronecker technique for improving quality of the signal in compressive sensing recovery

 This code shows the effect of Kronecker technique on compressive sensing recovery; this technique has been used for different signals and measurement matrices, 
  and it has been published in multiple journals and conference papers:

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
% For fast and efficient compression, sensing phase in compressive sensing can be done in very small size, because in this case: 
    %it requires very small measurement matrix, 
    %less number of multiplication and addition operations, and
    %less delay for generating compressed samples
    % which can be used for sensors with low computational resources.
 But sensing in small size degrades quality of recovery. Kronecker
 technique can be used in order to improve the quality of recovered signal.

%---------------------------------------CS Recovery algorithm
 Any recovery algorithm can be used. In this code, Sl0 which is very fast CS recovery algorithm is used.
 [sl0-reference]: http://ee.sharif.edu/~SLzero/

%----------------------------------Database
 In this code, an ECG signal from MIT-BIH Arrhythmia database which is public dataset is used.
 [Database Reference]:MIT-BIH Arrhythmia Database. [Online]. Available: http://www.physionet.org/physiobank/database/mitdb/
 

## Repo content
### Branches:
#### master
I will freely modify the code here to both test the pre-existing implementation (MATLAB), but also to port it to python for a project I'm working on.

Files:
- kronecker.ipynb (Jupyter Notebook):
     - python implementation

Directories:
- debugMatlab: 
     - modified version of original hadizand's code (MATLAB) with many debug prints to understand what happens at each step
     - added different files trying to understand and repeat the creation of wavelet dictionaries as performed in MATLAB through __wmpdictionary__

- hadizandMatlab:
     - original hadizand's code (MATLAB) without all the debug prints and added functions, just to clean test the hadizand implementation

#### hadizand-original-backup
Untouched copy of orginal folder


## How to use
### MATLAB content

You need to have MATLAB installed: `debugMATLAB`, `hadizandMATLAB` are stand-alone-content directories, just set MATLAB working directory to the one you want to run.
- Version: original code was written in 2008, I used it in 2024, no compatibility issues to this day. 

### Python
It slightly depends on the following choice you make: you can either choose to use only the provided RECORD 100.mat from the repository or use other records from the MIT-BIH Arrhythmia Database. 

#### Choice I: use only 100.mat in repo
- Today: 23 August 2024 (worked fine up to this date)
- Obviously need to have __python__ installed (__version 3.5 or above__, worked both with __3.8__ and __3.12__)
- You might want to install a python virtual environment, here shown with venv (but you're free to use _conda_ or whatever you prefer)
```sh
# virtual environment
python -m venv .kroneckerVenv
source .kroneckerVenv/bin/activate
```
- If working with Jupyter Notebook (.ipynb)
     - NOTE: to this day, I only provide a version in jupyter, my intention is to transpose to individual scripts once I have the time, for better portability and mantainability.
```sh
# ipykernel installation
pip install jupyter ipykernel
python -m ipykernel install --user --name=kroneckerVenv --display-name "Python (kroneckerVenv)"
```
- Install required python modules
```sh
# other modules
pip install -r .requirements.txt
```

#### Choice II: use freely [MIT-BIH Arrhythmia Database](https://physionet.org/content/mitdb/1.0.0/)
- There is a 