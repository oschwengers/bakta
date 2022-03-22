import json

from pathlib import Path
from subprocess import run

import bakta.config as cfg
import bakta.expert.protein_sequences as exp_aa_seq

import pytest


SEQUENCE = 'MRADEEPGDLSAVAQDYLKVIWTAQEWSQDKVSTKMLAERIGVSASTASESIRKLAEQGLVDHEKYGAVTLTDSGRRAALAMVRRHRLLETFLVNELGYRWDEVHDEA'


aa_min = {
    'id': 'min',
    'description': '~~~product~~~',
    'sequence': SEQUENCE
}
aa_min_gene = {
    'id': 'min',
    'description': 'gene~~~product~~~',
    'sequence': SEQUENCE
}
aa_min_dbxref = {
    'id': 'min',
    'description': 'gene~~~product~~~db-1:id-1',
    'sequence': SEQUENCE
}
aa_min_dbxrefs = {
    'id': 'min',
    'description': 'gene~~~product~~~db-1:id-1,db-2:id-2',
    'sequence': SEQUENCE
}
aa_full = {
    'id': 'full',
    'description': '90.0~~~80.0~~~80.0~~~gene~~~product~~~db-1:id-1,db-2:id-2',
    'sequence': SEQUENCE
}


aa_wrong_1 = {
    'id': 'low-cols',
    'description': '~~~product',
    'sequence': SEQUENCE
}
aa_wrong_2 = {
    'id': 'high-cols',
    'description': '90~~~80~~~80~~~gene~~~product~~~dbxref:dbxref~~~',
    'sequence': SEQUENCE
}
aa_wrong_3 = {
    'id': 'no-product',
    'description': 'gene~~~~~~dbxref:dbxref',
    'sequence': SEQUENCE
}
aa_wrong_4 = {
    'id': 'no-product-full',
    'description': '90~~~80~~~80~~~gene~~~~~~dbxref:dbxref',
    'sequence': SEQUENCE
}
aa_wrong_5 = {
    'id': 'wrong-dbxref',
    'description': 'gene~~~product~~~dbxrefdbxref',
    'sequence': SEQUENCE
}
aa_wrong_6 = {
    'id': 'wrong-dbxref-full',
    'description': '90~~~80~~~80~~~gene~~~product~~~dbxrefdbxref',
    'sequence': SEQUENCE
}
aa_wrong_7 = {
    'id': 'wrong-id',
    'description': 'ninety~~~80~~~80~~~gene~~~product~~~dbxref:dbxref',
    'sequence': SEQUENCE
}
aa_wrong_8 = {
    'id': 'wrong-min-query-cov',
    'description': '90~~~eighty~~~80~~~gene~~~product~~~dbxref:dbxref',
    'sequence': SEQUENCE
}
aa_wrong_9 = {
    'id': 'wrong-min-model-cov',
    'description': '90~~~80~~~eighty~~~gene~~~product~~~dbxref:dbxref',
    'sequence': SEQUENCE
}


@pytest.mark.parametrize(
    "parameters",
    [
        (aa_wrong_1),
        (aa_wrong_2),
        (aa_wrong_3),
        (aa_wrong_4),
        (aa_wrong_5),
        (aa_wrong_6),
        (aa_wrong_7),
        (aa_wrong_8),
        (aa_wrong_9)
    ]
)
def test_wrong_user_proteins_io(tmpdir, parameters):
    tmpdir = Path(tmpdir)
    cfg.user_proteins = tmpdir.joinpath('user.faa')
    write_tmp_faa(parameters, cfg.user_proteins)

    user_proteins_path = tmpdir.joinpath('user-clean.faa')
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        exp_aa_seq.write_user_protein_sequences(user_proteins_path)
    assert pytest_wrapped_e.type == SystemExit


@pytest.mark.parametrize(
    "parameters",
    [
        (aa_min),
        (aa_min_gene),
        (aa_min_dbxref),
        (aa_min_dbxrefs),
        (aa_full)
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


@pytest.mark.slow
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
            'bin/bakta', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', '--proteins', parameters,
            '--skip-tmrna', '--skip-trna', '--skip-rrna', '--skip-ncrna', '--skip-ncrna-region', '--skip-crispr', '--skip-sorf', '--skip-ori', '--skip-gap', 
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
        if('expert' in feat and 'user_proteins' in feat['expert']):
            user_prot_feats.append(feat)
        for db_xref in feat['db_xrefs']:
            if('EC' in db_xref):
                ec_annotated += 1
    assert ec_annotated == 1
    assert len(user_prot_feats) == 1
