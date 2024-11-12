
from os import path
from setuptools import setup, find_packages
import bakta


# Get the long description from the README file
setup_dir = path.abspath(path.dirname(__file__))
with open(path.join(setup_dir, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='bakta',
    version=bakta.__version__,
    description='Bakta: rapid & standardized annotation of bacterial genomes, MAGs & plasmids',
    keywords=['bioinformatics', 'annotation', 'bacteria', 'plasmids'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='GPLv3',
    author='Oliver Schwengers',
    author_email='oliver.schwengers@computational.bio.uni-giessen.de',
    url='https://github.com/oschwengers/bakta',
    packages=find_packages(include=['bakta', 'bakta.*']),
    python_requires='>=3.9, <3.12',
    include_package_data=False,
    zip_safe=False,
    install_requires=[
        'biopython >= 1.78',
        'xopen >= 1.5.0',
        'requests >= 2.25.1',
        'alive-progress >= 3.0.1',
        'PyYAML >= 6.0',
        'pyrodigal >= 3.5.0',
        'pyhmmer >= 0.10.15',
        'pycirclize >= 1.7.0'
    ],
    entry_points={
        'console_scripts': [
            'bakta=bakta.main:main',
            'bakta_proteins=bakta.proteins:main',
            'bakta_db=bakta.db:main',
            'bakta_plot=bakta.plot:main',
            'bakta_io=bakta.json_io:main'
        ]
    },
    classifiers=[
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English'
    ],
    project_urls={
        'Documentation': 'https://bakta.readthedocs.io',
        'Source': 'https://github.com/oschwengers/bakta',
        'Bug Reports': 'https://github.com/oschwengers/bakta/issues',
        'CI': 'https://github.com/oschwengers/bakta/actions'
    },
)
