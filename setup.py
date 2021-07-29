
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
    description='Bakta: rapid & comprehensive annotation of bacterial genomes & plasmids',
    keywords=['bioinformatics', 'annotation', 'bacteria', 'plasmids'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='GPLv3',
    author='Oliver Schwengers',
    author_email='oliver.schwengers@computational.bio.uni-giessen.de',
    url='https://github.com/oschwengers/bakta',
    packages=find_packages(include=['bakta', 'bakta.*']),
    python_requires='>=3.8',
    include_package_data=False,
    zip_safe=False,
    install_requires=[
        'biopython >= 1.78',
        'xopen >= 1.1.0',
        'requests >= 2.25.1',
        'alive-progress >= 1.6.2'
    ],
    entry_points={
        'console_scripts': [
            'bakta=bakta.main:main',
            'bakta_db=bakta.db:main'
        ]
    },
    classifiers=[
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.8',
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English'
    ],
    project_urls={
        'Bug Reports': 'https://github.com/oschwengers/bakta/issues',
        'Source': 'https://github.com/oschwengers/bakta'
    },
)
