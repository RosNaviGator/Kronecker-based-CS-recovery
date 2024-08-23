# CSV loggign branch (ipynb)
In this branch...
- Python project is all in a single Jupyter Notebook `kronecker.ipynb` 
- Both MATLAB and Python code are modified for __csv logging__: every important array is saved in a csv so that __you can compare the two implementations__ (at the end of python notebook there is code snippet that automatically performs the comparison)


## How to use (same as master ...)

---

### MATLAB content

You need to have MATLAB installed: `debugMATLAB`, `hadizandMATLAB` are stand-alone-content directories, just set MATLAB working directory to the one you want to run.
- Version: original code was written in 2008, I used it in 2024, no compatibility issues to this day. 

---

### Python
It slightly depends on the following choice you make: you can either choose to use only the provided RECORD 100.mat from the repository or use other records from the MIT-BIH Arrhythmia Database.

#### Python version and modules
- Today: 23 August 2024 (worked fine up to this date)
- Version: obviously need to have __python__ installed (__version 3.5 or above__, worked both with __3.8__ and __3.12__)

- Run all commands from `root` directory of the project

- You might want to install a python virtual environment, here shown with venv (but you're free to use _conda_ or whatever you prefer)
- Jupyter Notebook (.ipynb)
```sh
# virtual environment
python -m venv .kroneckerVenv
source .kroneckerVenv/bin/activate
```
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

#### Automatic installation script (alternative to previous commands)
```sh
source install.sh
```


# Kronecker-based-CS-recovery -- ECG data -- RosNaviGator
### Fork to study and to port to python (if I manage) the work of [Kronecker-based-CS-recovery](https://github.com/hadizand/Kronecker-based-CS-recovery)
