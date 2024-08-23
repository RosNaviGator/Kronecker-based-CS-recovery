# venv
python -m venv kroneckerVenv
source kroneckerVenv/bin/activate

# ipykernel installation
pip install jupyter ipykernel
python -m ipykernel install --user --name=kroneckerVenv --display-name "Python (kroneckerVenv)"

# other modules
pip install -r .requirements.txt
