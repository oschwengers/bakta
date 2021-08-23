import hashlib
import json
import logging
import os
import shutil
import subprocess as sp
import stat
import sys
import tarfile
import tempfile

from pathlib import Path

from alive_progress import alive_bar
import requests

import bakta
import bakta.constants as bc
import bakta.utils as bu


log = logging.getLogger('DB')


FILE_NAMES = [
        'antifam.h3f',
        'antifam.h3i',
        'antifam.h3m',
        'antifam.h3p',
        'bakta.db',
        'expert-protein-sequences.dmnd',
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
        'pfam.h3f',
        'pfam.h3i',
        'pfam.h3m',
        'pfam.h3p',
        'psc.dmnd',
        'rfam-go.tsv',
        'rRNA.i1f',
        'rRNA.i1i',
        'rRNA.i1m',
        'rRNA.i1p',
        'sorf.dmnd'
    ]


def check(db_path):
    """Check if database directory exists, is accessible and contains necessary files."""

    if(db_path is None):
        log.error('directory neither provided nor detected!')
        sys.exit('ERROR: database directory not provided nor detected! Please provide a valid path to the database directory.')

    if(not os.access(str(db_path), os.R_OK & os.X_OK)):
        log.error('directory (%s) not readable/accessible!', db_path)
        sys.exit(f'ERROR: database directory ({db_path}) not readable/accessible!')

    version_path = db_path.joinpath('version.json')
    if(not os.access(str(version_path), os.R_OK) or not version_path.is_file()):
        log.error('version file not readable!')
        sys.exit('ERROR: database version file (version.json) not readable!')

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

    log.info('detected: major=%i, minor=%i, date=%s', db_info['major'], db_info['minor'], db_info['date'])
    if(db_info['major'] < bakta.__db_schema_version__):
        log.error('wrong database version detected! required=%i, detected=%i', bakta.__db_schema_version__, db_info['major'])
        sys.exit(f"ERROR: wrong database version detected!\nBakta version {bakta.__version__} requires database version {bakta.__db_schema_version__}.x, but {db_info['major']}.{db_info['minor']} was detected. Please, update the database from https://doi.org/10.5281/zenodo.4247253")
    elif(db_info['major'] > bakta.__db_schema_version__):
        log.error('wrong database version detected! required=%i, detected=%i', bakta.__db_schema_version__, db_info['major'])
        sys.exit(f"ERROR: wrong database version detected!\nBakta version {bakta.__version__} requires database version {bakta.__db_schema_version__}.x, but {db_info['major']}.{db_info['minor']} was detected. Please, update Bakta or download a compatible database version from https://doi.org/10.5281/zenodo.4247253")

    for file_name in FILE_NAMES:
        path = db_path.joinpath(file_name)
        if(not os.access(str(path), os.R_OK) or not path.is_file()):
            log.error('file not readable! file=%s', file_name)
            sys.exit(f'ERROR: database file ({file_name}) not readable!')

    return db_info


def fetch_db_versions():
    try:
        with requests.get(bc.DB_VERSIONS_URL) as resp:
            versions = json.loads(resp.content)
    except IOError as e:
        print(e, file=sys.stderr)
        raise e
    else:
        return versions


def download(db_url, tarball_path):
    try:
        with tarball_path.open('wb') as fh_out, requests.get(db_url, stream=True) as resp:
            total_length = resp.headers.get('content-length')
            if(total_length is None):  # no content length header
                with alive_bar() as bar:
                    for data in resp.iter_content(chunk_size=1024*1024):
                        fh_out.write(data)
                        bar()
            else:
                total_length = int(int(total_length)/1024)  # length in Kb
                with alive_bar(total=total_length) as bar:
                    for data in resp.iter_content(chunk_size=1024*1024):
                        fh_out.write(data)
                        bar(incr=len(data)/1024)
    except IOError:
        sys.exit(f'ERROR: Could not download file from Zenodo! url={db_url}, path={tarball_path}')


def calc_md5_sum(tarball_path, buffer_size=1024*1024):
    md5 = hashlib.md5()
    with tarball_path.open('rb') as fh:
        data = fh.read(buffer_size)
        while data:
            md5.update(data)
            data = fh.read(buffer_size)
    return md5.hexdigest()


def untar(tarball_path, output_path):
    try:
        with tarball_path.open('rb') as fh_in, tarfile.open(fileobj=fh_in, mode='r:gz') as tar_file:
            tar_file.extractall(path=str(output_path))
    except OSError:
        sys.exit(f'ERROR: Could not extract {tarball_path} to {output_path}')


def main():
    # parse options and arguments
    parser = bu.init_parser()
    group_runtime = parser.add_argument_group('Runtime & auxiliary options')
    group_runtime.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    group_runtime.add_argument('--version', '-V', action='version', version=f'%(prog)s {bakta.__version__}')

    subparsers = parser.add_subparsers(dest='subcommand', help='sub-command help')
    parser_list = subparsers.add_parser('list', help='List available database versions')  # add list sub-command options
    parser_list.add_argument('--all', action='store_true', help='Show all versions including incompatible')

    parser_download = subparsers.add_parser('download', help='Download a database')  # add download sub-command options
    parser_download.add_argument('--output', '-o', action='store', default=Path.cwd(), help='output directory (default = current working directory)')
    parser_download.add_argument('--minor', '-n', action='store', type=int, default=0, help='Database minor version (default = most recent db minor version)')

    parser_update = subparsers.add_parser('update', help='Update an existing database to the most recent compatible version')  # add download sub-command options
    parser_update.add_argument('--db', '-d', action='store', default=None, help='Current database path (default = <bakta_path>/db). Can also be provided as BAKTA_DB environment variable.')
    parser_update.add_argument('--tmp-dir', '-t', action='store', dest='tmp_dir', default=Path.cwd(), help='Temporary directory to download & extract (default = current working directory)')

    args = parser.parse_args()
    if(args.subcommand == 'list'):
        print(f'Required database schema version: {bakta.__db_schema_version__}\n')
        versions = fetch_db_versions()
        if(not args.all):
            versions = [v for v in versions if v['major'] == bakta.__db_schema_version__]

        print('Available DB versions:')
        for v in sorted(versions, key=lambda v: (v['major'], v['minor'])):
            print(f"{v['major']}.{v['minor']}\t{v['date']}\t{v['doi']}")
    elif(args.subcommand == 'download'):
        bu.test_dependency(bu.DEPENDENCY_AMRFINDERPLUS)
        try:
            output_path = Path(args.output)
            if(not output_path.exists()):
                output_path.mkdir(parents=True, exist_ok=True)
            elif(not os.access(str(output_path), os.X_OK)):
                sys.exit(f'ERROR: output path ({output_path}) not accessible!')
            elif(not os.access(str(output_path), os.W_OK)):
                sys.exit(f'ERROR: output path ({output_path}) not writable!')
            output_path = output_path.resolve()
        except:
            sys.exit(f'ERROR: could not resolve or create output directory ({args.output})!')

        print('fetch DB versions...')
        versions = fetch_db_versions()
        compatible_versions = [v for v in versions if v['major'] == bakta.__db_schema_version__]
        if(len(compatible_versions) == 0):
            sys.exit(f'Error: no compatible version available for current major db version {bakta.__db_schema_version__}')
        else:
            print(f'\t... compatible DB versions: {len(compatible_versions)}')

        required_version = None
        if(args.minor > 0):
            for v in versions:
                if(v['minor'] == args.minor):
                    required_version = v
                    break
            if(required_version is None):
                sys.exit(f"requested DB minor version {args.minor} is not available. Please use 'bakta_db list' to get a list of available DB versions")
        else:
            compatible_sorted = sorted(compatible_versions, key=lambda v: v['minor'], reverse=True)
            required_version = compatible_sorted[0]

        tarball_path = output_path.joinpath('db.tar.gz')
        db_url = f"https://zenodo.org/record/{required_version['record']}/files/db.tar.gz"
        print(f"download database: v{required_version['major']}.{required_version['minor']}, {required_version['date']}, DOI: {required_version['doi']}, URL: {db_url}...")
        download(db_url, tarball_path)
        print('\t... done')

        print('check MD5 sum...')
        md5_sum = calc_md5_sum(tarball_path)
        if(md5_sum == required_version['md5']):
            print(f'\t...database file OK: {md5_sum}')
        else:
            sys.exit(f"Error: corrupt database file! MD5 should be '{required_version['md5']}' but is '{md5_sum}'")

        print(f'extract DB tarball: file={tarball_path}, output={output_path}')
        untar(tarball_path, output_path)
        tarball_path.unlink()

        db_path = output_path.joinpath('db')
        db_info = check(db_path)
        if(db_info['major'] != required_version['major']):
            sys.exit(f"ERROR: wrong major db detected! required={required_version['major']}, detected={db_info['major']}")
        elif(db_info['minor'] != required_version['minor']):
            sys.exit(f"ERROR: wrong minor db detected! required={required_version['minor']}, detected={db_info['minor']}")
        print('successfully downloaded Bakta database!')
        print(f"\tversion: {required_version['major']}.{required_version['minor']}")
        print(f"\tDOI: {required_version['doi']}")
        print(f'\tpath: {db_path}')

        try:
            db_path.chmod(stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)  # set write permissions on old (existing) directory with updated content
            for db_file_path in db_path.iterdir():
                db_file_path.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
        except:
            sys.exit(f'ERROR: cannot set read(|execute) permissions on new database! path={db_path}, owner={db_path.owner()}, group={db_path.group()}, permissions={oct(db_path.stat().st_mode )[-3:]}')

        print('update AMRFinderPlus database...')
        update_amrfinderplus_db(db_path)
        print('\t... done')

        print(f"\nRun Bakta using '--db {db_path}' or set a BAKTA_DB environment variable: 'export BAKTA_DB={db_path}'")
    elif(args.subcommand == 'update'):
        bu.test_dependency(bu.DEPENDENCY_AMRFINDERPLUS)
        env = os.environ.copy()
        if(args.db):
            db_dir = args.db
            try:
                db_tmp_path = Path(db_dir).resolve()
                if(db_tmp_path.is_dir()):
                    db_old_path = db_tmp_path
                    print(f'database provided via parameter: path={db_old_path}')
                else:
                    sys.exit(f'ERROR: unvalid database path! type=parameter, path={db_tmp_path}')
            except:
                sys.exit(f'ERROR: wrong database path! type=parameter, path={db_dir}')
        elif('BAKTA_DB' in env):
            db_dir = env['BAKTA_DB']
            try:
                db_tmp_path = Path(db_dir).resolve()
                if(db_tmp_path.is_dir()):
                    db_old_path = db_tmp_path
                    print(f'database provided via environment: path={db_old_path}')
                else:
                    sys.exit(f'ERROR: unvalid database path! type=environment, path={db_tmp_path}')
            except:
                sys.exit(f'ERROR: wrong database path! type=environment, BAKTA_DB={db_dir}')
        else:
            base_dir = Path(__file__).parent.parent
            db_tmp_path = base_dir.joinpath('db')
            if(db_tmp_path.is_dir()):
                db_old_path = db_tmp_path
                print(f'database detected in base-dir: path={db_old_path}')
            else:
                sys.exit('ERROR: database neither auto-detected nor provided!\nPlease, download the mandatory db and provide it either via the --db parameter, via a BAKTA_DB environment variable or copy it into the Bakta base directory.\nFor further information please read the readme.md')

        if(args.tmp_dir):
            tmp_path = Path(args.tmp_dir)
            if(not tmp_path.exists()):
                log.debug('dedicated temp dir does not exist! tmp-dir=%s', tmp_path)
                sys.exit(f'ERROR: dedicated temp dir ({tmp_path}) does not exist!')
            else:
                log.info('use dedicated temp dir: path=%s', tmp_path)
                tmp_path = Path(tempfile.mkdtemp(dir=str(tmp_path)))
        else:
            tmp_path = Path(tempfile.mkdtemp())

        db_old_info = check(db_old_path)
        print(f"existing database: v{db_old_info['major']}.{db_old_info['minor']}")
        print('fetch DB versions...')
        versions = fetch_db_versions()
        compatible_versions = [v for v in versions if v['major'] == bakta.__db_schema_version__]
        if(len(compatible_versions) == 0):
            sys.exit(f'Error: no compatible version available for current major db version {bakta.__db_schema_version__}')
        else:
            print(f'\t... compatible DB versions: {len(compatible_versions)}')

        compatible_sorted = sorted(compatible_versions, key=lambda v: v['minor'], reverse=True)
        if(compatible_sorted[0]['minor'] <= db_old_info['minor']):
            print(f"Database version {db_old_info['major']}.{db_old_info['minor']} is up-to-date")
            sys.exit()
        required_version = compatible_sorted[0]

        tarball_path = tmp_path.joinpath('db.tar.gz')
        db_url = f"https://zenodo.org/record/{required_version['record']}/files/db.tar.gz"
        print(f"download database: v{required_version['major']}.{required_version['minor']}, {required_version['date']}, DOI: {required_version['doi']}, URL: {db_url}...")
        download(db_url, tarball_path)
        print('\t... done')

        print('check MD5 sum...')
        md5_sum = calc_md5_sum(tarball_path)
        if(md5_sum == required_version['md5']):
            print(f'\t...database file OK: {md5_sum}')
        else:
            sys.exit(f"Error: corrupt database file! MD5 should be '{required_version['md5']}' but is '{md5_sum}'")

        print(f'extract DB tarball: file={tarball_path}, output-directory={tmp_path}')
        untar(tarball_path, tmp_path)
        tarball_path.unlink()

        db_new_path = tmp_path.joinpath('db')
        db_new_info = check(db_new_path)
        if(db_new_info['major'] != required_version['major']):
            sys.exit(f"ERROR: wrong major db detected! required={required_version['major']}, detected={db_new_info['major']}")
        elif(db_new_info['minor'] != required_version['minor']):
            sys.exit(f"ERROR: wrong minor db detected! required={required_version['minor']}, detected={db_new_info['minor']}")
        print('successfully downloaded Bakta DB:')
        print(f"\tversion: {required_version['major']}.{required_version['minor']}")
        print(f"\tDOI: {required_version['doi']}")
        print(f'\tpath: {db_new_path}')
        print('remove old database...')
        try:
            db_old_path.chmod(stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)  # set write permissions on old directory
            for db_old_file_path in db_old_path.iterdir():
                db_old_file_path.chmod(stat.S_IRUSR | stat.S_IWUSR)
        except:
            sys.exit(f'ERROR: cannot set read|write|execute permissions on old database! path={db_old_path}, owner={db_old_path.owner()}, group={db_old_path.group()}, permissions={oct(db_old_path.stat().st_mode )[-3:]}')
        try:
            shutil.rmtree(db_old_path)
        except:
            sys.exit(f'ERROR: cannot remove old database! path={db_old_path}, owner={db_old_path.owner()}, group={db_old_path.group()}, permissions={oct(db_old_path.stat().st_mode )[-3:]}')
        db_old_path.mkdir()

        try:
            db_new_path.chmod(stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)  # set write permissions on db_new_path directory
            for db_new_file_path in db_new_path.iterdir():
                db_new_file_path.chmod(stat.S_IRUSR | stat.S_IWUSR)
        except:
            sys.exit(f'ERROR: cannot set read|write|execute permissions on new database! path={db_new_path}, owner={db_new_path.owner()}, group={db_new_path.group()}, permissions={oct(db_new_path.stat().st_mode )[-3:]}')
        try:
            for db_new_file_path in db_new_path.iterdir():  # move new db files into old (existing) db directory
                file_name = db_new_file_path.name
                shutil.move(db_new_file_path, db_old_path.joinpath(file_name))
        except:
            sys.exit(f'ERROR: cannot move new database to existing path! new-path={db_new_path}, existing-path={db_old_path.parent}')
        shutil.rmtree(tmp_path)

        try:
            db_old_path.chmod(stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)  # set write permissions on old (existing) directory with updated content
            for db_old_file_path in db_old_path.iterdir():
                db_old_file_path.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
        except:
            sys.exit(f'ERROR: cannot set read(|execute) permissions on new database! path={db_old_path}, owner={db_old_path.owner()}, group={db_old_path.group()}, permissions={oct(db_old_path.stat().st_mode )[-3:]}')

        print('\t... done')

        print('update AMRFinderPlus database...')
        update_amrfinderplus_db(db_old_path)
        print('\t... done')

        print(f"\nRun Bakta using '--db {db_old_path}' or set a BAKTA_DB environment variable: 'export BAKTA_DB={db_old_path}'")
    else:
        parser.print_help()
        sys.exit('Error: no subcommand provided!')


def update_amrfinderplus_db(db_path):
    amrfinderplus_db_path = db_path.joinpath('amrfinderplus-db')
    cmd = [
        'amrfinder_update',
        '--database', str(amrfinderplus_db_path),
        '--force_update'
    ]
    log.debug('cmd=%s', cmd)
    proc = sp.run(
        cmd,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if(proc.returncode != 0):
        log.debug('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        log.warning('AMRFinderPlus failed! amrfinder-error-code=%d', proc.returncode)
        sys.exit(f"ERROR: AMRFinderPlus failed! command: 'amrfinder_update --force_update --database {amrfinderplus_db_path}', error code: {proc.returncode}")


if __name__ == '__main__':
    main()
