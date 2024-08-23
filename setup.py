from setuptools import setup, find_packages

setup(
    name='compSensPack',
    version='1.0.0',
    description='A library for compressed sensing with Kronecker Technique.',
    author='RosNaviGator',
    author_email='borndescalo@gmail.com',
    packages=find_packages(),  # Automatically find all packages and subpackages
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'PyWavelets',
    ],
    entry_points={
        'console_scripts': [
            'run-main-analysis = scripts.run_main_analysis:main',  # Adds a command line script
        ],
    },
)
