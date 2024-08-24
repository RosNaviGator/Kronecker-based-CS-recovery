# Kronecker-based-CS-recovery -- ECG data -- RosNaviGator


## Fork to study and to port to python (if I manage) the work of [Kronecker-based-CS-recovery](https://github.com/hadizand/Kronecker-based-CS-recovery)

I will leave as description the original one by Hadi Zanddizari, it also __includes all the references__ to __Kronecker Technique__, __Compressed Sensing__, __SL0 recovery algorithm__, __MIT-BIH Arrhythmia database__ 


### Original description (by [Hadi Zanddizari](https://github.com/hadizand))

**Kronecker Technique for Improving Signal Quality in Compressive Sensing Recovery**

This code demonstrates the effect of the Kronecker technique on compressive sensing recovery. This technique has been used for different signals and measurement matrices and has been published in multiple journals and conference papers:

1. **H. Zanddizari, S. Rajan, and H. Zarrabi**, “Increasing the quality of reconstructed signal in compressive sensing utilizing Kronecker technique,”  
   *Biomedical Engineering Letters*, vol. 8, no. 2, pp. 239–247, May 2018.

2. **D. Mitra, H. Zanddizari, and S. Rajan**, "Investigation of Kronecker-based recovery of compressed ECG signal,"  
   *IEEE Transactions on Instrumentation and Measurement*, pp. 1-1, 2019.

3. **D. Mitra, H. Zanddizari, and S. Rajan**, “Improvement of signal quality during recovery of compressively sensed ECG signals,”  
   *2018 IEEE International Symposium on Medical Measurements and Applications (MeMeA)*, June 2018, pp. 1–5.

4. **D. Mitra, H. Zanddizari, and S. Rajan**, “Improvement of recovery in segmentation-based parallel compressive sensing,”  
   *2018 IEEE International Symposium on Signal Processing and Information Technology (ISSPIT)*, Dec 2018, pp. 501–506.

**Original Author (MATLAB)**: Hadi Zanddizari  
Email: [hadiz@mail.usf.edu](mailto:hadiz@mail.usf.edu), [hadizand@alumni.iust.ac.ir](mailto:hadizand@alumni.iust.ac.ir)

#### The Main Objective of This Approach

For fast and efficient compression, the sensing phase in compressive sensing can be done in a very small size, because in this case:

- It requires a very small measurement matrix,
- Less number of multiplication and addition operations,
- Less delay for generating compressed samples,

These characteristics are particularly useful for sensors with low computational resources. However, sensing in small sizes degrades the quality of recovery. The Kronecker technique can be used to improve the quality of the recovered signal.

#### CS Recovery Algorithm

Any recovery algorithm can be used. In this code, **Sl0**, which is a very fast CS recovery algorithm, is used.  
SL0 Reference: [http://ee.sharif.edu/~SLzero/](http://ee.sharif.edu/~SLzero/)

#### Database

In this code, an ECG signal from the **MIT-BIH Arrhythmia Database**, which is a public dataset, is used.  
Database Reference: MIT-BIH Arrhythmia Database. Available: [http://www.physionet.org/physiobank/database/mitdb/](http://www.physionet.org/physiobank/database/mitdb/)

## Repo content

### Branches:

---

#### master
- Python implementation that __partially__ correctly emulates the original MATLAB code
- I never managed to build a __dwt dictionary__ in the way MATLAB's `wmpdictionary`, note that in the future Mathworks plans to remove the obsolete `wmpdictionary`, it will be substituted by `sensingDictionary`

---

#### csv-logging-ipynb
- Python code is all in a Jupyter Notebook
- Both `MATLAB` and `Python` codes save all relevant variables in `csv` so that the two languages implementation can be compared.

---

#### hadizand-original-backup
Untouched copy of orginal folder

---



## How to use

### MATLAB content

You need to have MATLAB installed: `debugMATLAB`, `hadizandMATLAB` are stand-alone-content directories, just set MATLAB working directory to the one you want to run.
- Version: original code was written in 2008, I used it in 2024, no compatibility issues to this day.
- `wmpdictionary` is obsolete, will be soon replaced by `sensingDictionary`

---

### Python
It slightly depends on the following choice you make: you can either choose to use only the provided RECORD 100.mat from the repository or use other records from the MIT-BIH Arrhythmia Database.

#### Python version and modules
- Today: 23 August 2024 (worked fine up to this date)
- Version: obviously need to have __python__ installed (__version 3.5 or above__, worked both with __3.8__ and __3.12__)

- Run all commands from `root` directory of the project

- You might want to install a python virtual environment, here shown with venv (but you're free to use _conda_ or whatever you prefer)
```sh
# virtual environment
python -m venv .kroneckerVenv
source .kroneckerVenv/bin/activate
```
- Install this project package `compSensPack` and all it's requirements
```sh
# install compSensPack
pip install .
```

#### Automatic installation script (alternative to previous ones)
```sh
source install.sh
```

#### Use freely whole [MIT-BIH Arrhythmia Database](https://physionet.org/content/mitdb/1.0.0/)
- There is a dedicated python module to use the MIT-BIH Arrhythmia Database
- __Beware:__ when I used this library in 2024 it __didn't work__ with latest version of python (__3.12__), it did work fine using __python3.8__
```sh
# Python module to use MIT-BIH databases bind to Physionet.org
pip install wfdb
```
- __I didn't use it__ for my implementation, RECORD 100 already present in original repo was enough for the purpose of my study, so _the rest is up to you..._