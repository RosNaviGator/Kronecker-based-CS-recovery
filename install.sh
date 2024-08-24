# venv
python -m venv kroneckerVenv
source kroneckerVenv/bin/activate

# ipykernel installation
pip install jupyter ipykernel
python -m ipykernel install --user --name=kroneckerVenv --display-name "Python (kroneckerVenv)"

# install compSensPack + all dependencies
pip install .