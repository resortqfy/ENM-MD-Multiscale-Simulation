from setuptools import setup, find_packages

setup(
    name='enm_md_coupling',
    version='0.1.0',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'numpy>=1.21.0',
        'scipy>=1.7.0',
        'matplotlib>=3.5.0',
        'pandas>=1.3.0',
        'scikit-learn>=1.0.0',
        'numba>=0.56.0',
        'MDAnalysis>=2.0.0',
        'biopython>=1.79',
        'seaborn>=0.11.0',
    ],
    tests_require=['pytest>=6.0.0'],
    description='ENM-MD Multiscale Protein Dynamics Simulation Framework',
    author='Your Name',
    author_email='your.email@example.com',
    license='MIT',
)

