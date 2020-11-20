
import json
import logging
import os
import sys

import bakta
import bakta.config as cfg

log = logging.getLogger('DB')


def check():
    """Check if database directory exists, is accessible and contains necessary files."""

    if(cfg.db_path is None):
        log.error('directory not provided nor detected!')
        sys.exit('ERROR: database directory not provided nor detected! Please provide a valid path to the database directory.')

    if(not os.access(str(cfg.db_path), os.R_OK & os.X_OK)):
        log.error('directory (%s) not readable/accessible!', cfg.db_path)
        sys.exit(f'ERROR: database directory ({cfg.db_path}) not readable/accessible!')

    version_path = cfg.db_path.joinpath('version.json')
    if(not os.access(str(version_path), os.R_OK) or not version_path.is_file()):
        log.error('version file not readable!')
        sys.exit('ERROR: database version file (db.yaml) not readable!')
    
    try:
        with version_path.open() as fh:
            db_info = json.load(fh)
    except:
        log.exception('could not parse database version file!')
        sys.exit('ERROR: could not parse database version file!')
    
    for key in ['date', 'major', 'minor']:
        if(key not in db_info):
            log.error('wrong db version info file content! missed key=%s', key)
            sys.exit(f"ERROR: wrong db version info file format! Missed key '{key}' in JSON structure.")
    
    log.info('detected: major=%i', db_info['major'])
    log.info('detected: minor=%i', db_info['minor'])
    log.info('detected: date=%s', db_info['date'])
    if(db_info['major'] < bakta.__db_schema_version__):
        log.error('wrong database version detected! required=%i, detected=%i', bakta.__db_schema_version__, db_info['major'])
        sys.exit(f"ERROR: wrong database version detected!\nBakta version {bakta.__version__} requires database version {bakta.__db_schema_version__}.x, but {db_info['major']}.{db_info['minor']} was detected. Please, update the database from https://doi.org/10.5281/zenodo.4247253")
    elif(db_info['major'] > bakta.__db_schema_version__):
        log.error('wrong database version detected! required=%i, detected=%i', bakta.__db_schema_version__, db_info['major'])
        sys.exit(f"ERROR: wrong database version detected!\nBakta version {bakta.__version__} requires database version {bakta.__db_schema_version__}.x, but {db_info['major']}.{db_info['minor']} was detected. Please, update Bakta or download a compatible database version from https://doi.org/10.5281/zenodo.4247253")

    file_names = [
        'antifam.h3f',
        'antifam.h3i',
        'antifam.h3m',
        'antifam.h3p',
        'bakta.db',
        'ncRNA-genes.i1f',
        'ncRNA-genes.i1i',
        'ncRNA-genes.i1m',
        'ncRNA-genes.i1p',
        'ncRNA-regions.i1f',
        'ncRNA-regions.i1i',
        'ncRNA-regions.i1m',
        'ncRNA-regions.i1p',
        'oric.fna',
        'orit.fna',
        'psc.dmnd',
        'rfam-go.tsv',
        'rRNA.i1f',
        'rRNA.i1i',
        'rRNA.i1m',
        'rRNA.i1p',
        'sorf.dmnd'
    ]

    for file_name in file_names:
        path = cfg.db_path.joinpath(file_name)
        if(not os.access(str(path), os.R_OK) or not path.is_file()):
            log.error('file not readable! file=%s', file_name)
            sys.exit(f'ERROR: database file ({file_name}) not readable!')

    return db_info
