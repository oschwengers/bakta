
from bakta import utils as bu


#  bu.check_version(tool_version, semver_depth, min_version, max_version, tool_name) == True/False


def test_major():
    # meeting conditions without max version
    assert bu.check_version(bu.Version(1,0,0), bu.Version(1,0,0), bu.Version(bu.Max_version,bu.Max_version,bu.Max_version), 'foo') == True
    assert bu.check_version(bu.Version(1,0), bu.Version(1,0,0), bu.Version(bu.Max_version,bu.Max_version,bu.Max_version), 'foo') == True
    assert bu.check_version(bu.Version(1), bu.Version(1,0,0), bu.Version(bu.Max_version,bu.Max_version,bu.Max_version), 'foo') == True

    # meeting conditions without min version
    assert bu.check_version(bu.Version(1,0,0), bu.Version(bu.Min_version,bu.Min_version,bu.Min_version), bu.Version(1,0,0), 'foo') == True
    assert bu.check_version(bu.Version(1,0), bu.Version(bu.Min_version,bu.Min_version,bu.Min_version), bu.Version(1,0,0), 'foo') == True
    assert bu.check_version(bu.Version(1), bu.Version(bu.Min_version,bu.Min_version,bu.Min_version), bu.Version(1,0,0), 'foo') == True

    # meeting conditions on lower bound
    assert bu.check_version(bu.Version(1,0,0), bu.Version(1,0,0), bu.Version(1,0,0), 'foo') == True
    assert bu.check_version(bu.Version(1,0), bu.Version(1,0,0), bu.Version(1,0), 'foo') == True
    assert bu.check_version(bu.Version(1), bu.Version(1,0,0), bu.Version(1), 'foo') == True
    
    # meeting conditions on upper bound
    assert bu.check_version(bu.Version(2,2,2), bu.Version(1,0,0), bu.Version(2,2,2), 'foo') == True
    assert bu.check_version(bu.Version(2,2), bu.Version(1,1), bu.Version(2,2), 'foo') == True
    assert bu.check_version(bu.Version(2), bu.Version(1), bu.Version(2), 'foo') == True
    
    # meeting conditions in between
    assert bu.check_version(bu.Version(1,0,1), bu.Version(1,0,0), bu.Version(2,0,0), 'foo') == True
    assert bu.check_version(bu.Version(1,0,100), bu.Version(1,0,0), bu.Version(2,0,0), 'foo') == True
    assert bu.check_version(bu.Version(1,1,0), bu.Version(1,0,0), bu.Version(2,0,0), 'foo') == True
    assert bu.check_version(bu.Version(1,100,0), bu.Version(1,0,0), bu.Version(2,0,0), 'foo') == True

    assert bu.check_version(bu.Version(1,1), bu.Version(1,0), bu.Version(2,0), 'foo') == True
    assert bu.check_version(bu.Version(1,100), bu.Version(1,0), bu.Version(2,0), 'foo') == True
    assert bu.check_version(bu.Version(1,1), bu.Version(1,0), bu.Version(2,0), 'foo') == True
    assert bu.check_version(bu.Version(1,100), bu.Version(1,0), bu.Version(2,0), 'foo') == True
    
    # meeting conditions in between
    assert bu.check_version(bu.Version(2,1,1), bu.Version(1,2,10), bu.Version(3,0,0), 'foo') == True
    assert bu.check_version(bu.Version(2,2,1), bu.Version(1,2,10), bu.Version(3,0,0), 'foo') == True
    
    
    # unmeeting major conditions without max version
    assert bu.check_version(bu.Version(1,0,0), bu.Version(2,0,0), bu.Version(bu.Max_version,bu.Max_version,bu.Max_version), 'foo') == False
    assert bu.check_version(bu.Version(1,0), bu.Version(2,0), bu.Version(bu.Max_version,bu.Max_version,bu.Max_version), 'foo') == False
    assert bu.check_version(bu.Version(1), bu.Version(2), bu.Version(bu.Max_version,bu.Max_version,bu.Max_version), 'foo') == False

    # unmeeting major conditions without min version
    assert bu.check_version(bu.Version(2,0,0), bu.Version(bu.Min_version,bu.Min_version,bu.Min_version), bu.Version(1,0,0), 'foo') == False
    assert bu.check_version(bu.Version(2,0), bu.Version(bu.Min_version,bu.Min_version,bu.Min_version), bu.Version(1,0,0), 'foo') == False
    assert bu.check_version(bu.Version(2), bu.Version(bu.Min_version,bu.Min_version,bu.Min_version), bu.Version(1,0,0), 'foo') == False

    # unmeeting major conditions on lower bound
    assert bu.check_version(bu.Version(1,0,0), bu.Version(2,0,0), bu.Version(3,0,0), 'foo') == False
    assert bu.check_version(bu.Version(1,0), bu.Version(2,0), bu.Version(3,0), 'foo') == False
    assert bu.check_version(bu.Version(1), bu.Version(2), bu.Version(3), 'foo') == False
    
    # unmeeting major conditions on upper bound
    assert bu.check_version(bu.Version(3,0,0), bu.Version(1,0,0), bu.Version(2,0,0), 'foo') == False
    assert bu.check_version(bu.Version(3,0), bu.Version(1,1), bu.Version(2,0), 'foo') == False
    assert bu.check_version(bu.Version(3), bu.Version(1), bu.Version(2), 'foo') == False


def test_minor():
    # meeting conditions without max version
    assert bu.check_version(bu.Version(1,10,0), bu.Version(1,1,0), bu.Version(bu.Max_version,bu.Max_version,bu.Max_version), 'foo') == True
    assert bu.check_version(bu.Version(1,10), bu.Version(1,1,0), bu.Version(bu.Max_version,bu.Max_version,bu.Max_version), 'foo') == True

    # meeting conditions without min version
    assert bu.check_version(bu.Version(1,1,0), bu.Version(bu.Min_version,bu.Min_version,bu.Min_version), bu.Version(1,10,0), 'foo') == True
    assert bu.check_version(bu.Version(1,1), bu.Version(bu.Min_version,bu.Min_version,bu.Min_version), bu.Version(1,10,0), 'foo') == True

    # meeting conditions on lower bound
    assert bu.check_version(bu.Version(1,0,0), bu.Version(1,0,0), bu.Version(1,10,0), 'foo') == True
    assert bu.check_version(bu.Version(1,0), bu.Version(1,0), bu.Version(1,10), 'foo') == True
    
    # meeting conditions on upper bound
    assert bu.check_version(bu.Version(1,10,0), bu.Version(1,0,0), bu.Version(1,10,0), 'foo') == True
    assert bu.check_version(bu.Version(1,10), bu.Version(1,0), bu.Version(1,10), 'foo') == True
    
    # meeting conditions in between
    assert bu.check_version(bu.Version(1,5,0), bu.Version(1,0,0), bu.Version(1,10,0), 'foo') == True
    

    # unmeeting major conditions without max version
    assert bu.check_version(bu.Version(1,0,0), bu.Version(1,10,0), bu.Version(bu.Max_version,bu.Max_version,bu.Max_version), 'foo') == False
    assert bu.check_version(bu.Version(1,0), bu.Version(1,10), bu.Version(bu.Max_version,bu.Max_version,bu.Max_version), 'foo') == False

    # unmeeting major conditions without min version
    assert bu.check_version(bu.Version(1,10,0), bu.Version(bu.Min_version,bu.Min_version,bu.Min_version), bu.Version(1,0,0), 'foo') == False
    assert bu.check_version(bu.Version(1,10), bu.Version(bu.Min_version,bu.Min_version,bu.Min_version), bu.Version(1,0), 'foo') == False

    # unmeeting major conditions on lower bound
    assert bu.check_version(bu.Version(1,0,0), bu.Version(1,10,0), bu.Version(1,0,0), 'foo') == False
    assert bu.check_version(bu.Version(1,0), bu.Version(1,10), bu.Version(1,0), 'foo') == False
    
    # unmeeting major conditions on upper bound
    assert bu.check_version(bu.Version(1,10,0), bu.Version(1,0,0), bu.Version(1,5,0), 'foo') == False
    assert bu.check_version(bu.Version(1,10), bu.Version(1,0), bu.Version(1,5), 'foo') == False


def test_patch():
    # meeting conditions without max version
    assert bu.check_version(bu.Version(1,0,10), bu.Version(1,0,1), bu.Version(bu.Max_version,bu.Max_version,bu.Max_version), 'foo') == True

    # meeting conditions without min version
    assert bu.check_version(bu.Version(1,0,10), bu.Version(bu.Min_version,bu.Min_version,bu.Min_version), bu.Version(1,0,10), 'foo') == True

    # meeting conditions on lower bound
    assert bu.check_version(bu.Version(1,0,1), bu.Version(1,0,1), bu.Version(1,0,10), 'foo') == True
    
    # meeting conditions on upper bound
    assert bu.check_version(bu.Version(1,0,10), bu.Version(1,0,0), bu.Version(1,0,10), 'foo') == True
    
    # meeting conditions in between
    assert bu.check_version(bu.Version(1,0,5), bu.Version(1,0,1), bu.Version(1,0,10), 'foo') == True
    

    # unmeeting major conditions without max version
    assert bu.check_version(bu.Version(1,0,0), bu.Version(1,0,10), bu.Version(bu.Max_version,bu.Max_version,bu.Max_version), 'foo') == False

    # unmeeting major conditions without min version
    assert bu.check_version(bu.Version(1,0,10), bu.Version(bu.Min_version,bu.Min_version,bu.Min_version), bu.Version(1,0,0), 'foo') == False

    # unmeeting major conditions on lower bound
    assert bu.check_version(bu.Version(1,0,0), bu.Version(1,0,10), bu.Version(1,0,10), 'foo') == False
    
    # unmeeting major conditions on upper bound
    assert bu.check_version(bu.Version(1,0,10), bu.Version(1,0,0), bu.Version(1,0,5), 'foo') == False
