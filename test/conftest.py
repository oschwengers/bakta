from _pytest.mark import Mark


global FILES, SKIP_PARAMETERS

FILES = [
    'test.log',
    'test.json',
    'test.tsv',
    'test.gff3',
    'test.gbff',
    'test.embl',
    'test.fna',
    'test.faa'
]

SKIP_PARAMETERS = [
    '--skip-tmrna',
    '--skip-trna',
    '--skip-rrna',
    '--skip-ncrna',
    '--skip-ncrna-region',
    '--skip-crispr',
    '--skip-cds',
    '--skip-sorf',
    '--skip-ori',
    '--skip-gap'
]


empty_mark = Mark('', [], {})


def by_slow_marker(item):
    return item.get_closest_marker('slow', default=empty_mark)


def pytest_collection_modifyitems(items):
    items.sort(key=by_slow_marker, reverse=False)
