from subprocess import run


def setup_package():
    """Update AMRFinderPlus db."""
    proc = run(['amrfinder', '--update'])
    assert proc.returncode == 0