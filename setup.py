from setuptools import setup, find_packages

VERSION = '0.99.4a0'
URL = 'https://github.com/galaxyproteomics/metaquantome'
AUTHOR = 'Caleb Easterly'
AUTHOR_EMAIL = 'caleb.easterly@gmail.com'
setup(
    name='metaquantome',
    version=VERSION,
    packages=find_packages(),
    package_data={'metaquantome': ['data/*', 'modules/viz.R']},
    include_package_data=True,
    license='Apache License 2.0',
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    description='Quantitative metaproteomics analysis',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research'
    ],
    install_requires=[
        'pandas',
        'ete3',
        'goatools',
        'numpy',
        'statsmodels',
        'biopython'
    ],
    python_requires='>=3.5',
    entry_points={
        'console_scripts': ['metaquantome=metaquantome.cli:cli'],
    },
    project_urls={
        "Bug Tracker": "https://github.com/galaxyproteomics/metaquantome/issues",
        "Source Code": "https://github.com/galaxyproteomics/metaquantome",
    }
)