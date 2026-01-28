import json

from pathlib import Path
from subprocess import run

import bakta.config as cfg
import bakta.expert.protein_sequences as exp_aa_seq

import pytest


SEQUENCE = 'MRADEEPGDLSAVAQDYLKVIWTAQEWSQDKVSTKMLAERIGVSASTASESIRKLAEQGLVDHEKYGAVTLTDSGRRAALAMVRRHRLLETFLVNELGYRWDEVHDEA'


AA_MIN = {
    'id': 'min',
    'description': '~~~product~~~',
    'sequence': SEQUENCE
}
AA_MIN_GENE = {
    'id': 'min',
    'description': 'gene~~~product~~~',
    'sequence': SEQUENCE
}
AA_MIN_DBXREF = {
    'id': 'min',
    'description': 'gene~~~product~~~db-1:id-1',
    'sequence': SEQUENCE
}
AA_MIN_DBXREFS = {
    'id': 'min',
    'description': 'gene~~~product~~~db-1:id-1,db-2:id-2',
    'sequence': SEQUENCE
}
AA_FULL = {
    'id': 'full',
    'description': '90.0~~~80.0~~~80.0~~~gene~~~product~~~db-1:id-1,db-2:id-2',
    'sequence': SEQUENCE
}


@pytest.mark.parametrize(
    "parameters",
    [
        AA_MIN,
        AA_MIN_GENE,
        AA_MIN_DBXREF,
        AA_MIN_DBXREFS,
        AA_FULL
    ]
)
def test_user_proteins_io(parameters, tmpdir):
    tmpdir = Path(tmpdir)
    cfg.user_proteins = tmpdir.joinpath('user.faa')
    write_tmp_faa(parameters, cfg.user_proteins)

    user_proteins_path = tmpdir.joinpath('user-clean.faa')
    exp_aa_seq.write_user_protein_sequences(user_proteins_path)


def write_tmp_faa(aa, aa_path):
    with aa_path.open('w') as fh:
        fh.write(f">{aa['id']} {aa['description']}\n{aa['sequence']}\n")


@pytest.mark.parametrize(
    "parameters",
    [
        'test/data/user-proteins.faa',
        'test/data/user-proteins.faa.gz',
        'test/data/user-proteins.gbff',
        'test/data/user-proteins.gbff.gz'
    ]
)
def test_user_proteins(parameters, tmpdir):
    # fast test skipping all feature detections
    proc = run(
        [
            'bin/bakta', '--db', 'test/db', '--output', tmpdir, '--force', '--prefix', 'test', '--proteins', parameters,
            '--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-sorf', '--skip-ori', '--skip-gap', '--skip-plot',
            'test/data/NC_002127.1.fna'
        ]
    )
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    results_path = Path(tmpdir_path.joinpath('test.json'))
    assert Path.exists(results_path)
    results = None
    with results_path.open() as fh:
        results = json.load(fh)
    assert results is not None
    
    ec_annotated = 0
    user_prot_feats = []
    for feat in results['features']:
        for expert in feat.get('expert', []):
            if(expert['type'] == 'user_proteins'):
                user_prot_feats.append(feat)
                ec_annotated += len([x for x in feat['db_xrefs'] if 'EC' in x])
    assert ec_annotated == 1
    assert len(user_prot_feats) == 1
