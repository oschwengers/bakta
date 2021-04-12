
from bakta import utils as bu


#  bu.check_version(tool_version, semver_depth, VERSION_MIN_DIGIT, VERSION_MAX_DIGIT, tool_name) == True/False


def test_major():
    # meeting conditions without max version
    assert bu.check_version(bu.Version(1,0,0), bu.Version(1,0,0), bu.Version(bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT)) == True
    assert bu.check_version(bu.Version(1,0), bu.Version(1,0,0), bu.Version(bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT)) == True
    assert bu.check_version(bu.Version(1), bu.Version(1,0,0), bu.Version(bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT)) == True

    # meeting conditions without min version
    assert bu.check_version(bu.Version(1,0,0), bu.Version(bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT), bu.Version(1,0,0)) == True
    assert bu.check_version(bu.Version(1,0), bu.Version(bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT), bu.Version(1,0,0)) == True
    assert bu.check_version(bu.Version(1), bu.Version(bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT), bu.Version(1,0,0)) == True

    # meeting conditions on lower bound
    assert bu.check_version(bu.Version(1,0,0), bu.Version(1,0,0), bu.Version(1,0,0)) == True
    assert bu.check_version(bu.Version(1,0), bu.Version(1,0,0), bu.Version(1,0)) == True
    assert bu.check_version(bu.Version(1), bu.Version(1,0,0), bu.Version(1)) == True
    
    # meeting conditions on upper bound
    assert bu.check_version(bu.Version(2,2,2), bu.Version(1,0,0), bu.Version(2,2,2)) == True
    assert bu.check_version(bu.Version(2,2), bu.Version(1,1), bu.Version(2,2)) == True
    assert bu.check_version(bu.Version(2), bu.Version(1), bu.Version(2)) == True
    
    # meeting conditions in between
    assert bu.check_version(bu.Version(1,0,1), bu.Version(1,0,0), bu.Version(2,0,0)) == True
    assert bu.check_version(bu.Version(1,0,100), bu.Version(1,0,0), bu.Version(2,0,0)) == True
    assert bu.check_version(bu.Version(1,1,0), bu.Version(1,0,0), bu.Version(2,0,0)) == True
    assert bu.check_version(bu.Version(1,100,0), bu.Version(1,0,0), bu.Version(2,0,0)) == True

    assert bu.check_version(bu.Version(1,1), bu.Version(1,0), bu.Version(2,0)) == True
    assert bu.check_version(bu.Version(1,100), bu.Version(1,0), bu.Version(2,0)) == True
    assert bu.check_version(bu.Version(1,1), bu.Version(1,0), bu.Version(2,0)) == True
    assert bu.check_version(bu.Version(1,100), bu.Version(1,0), bu.Version(2,0)) == True
    
    # meeting conditions in between
    assert bu.check_version(bu.Version(2,1,1), bu.Version(1,2,10), bu.Version(3,0,0)) == True
    assert bu.check_version(bu.Version(2,2,1), bu.Version(1,2,10), bu.Version(3,0,0)) == True
    
    
    # unmeeting major conditions without max version
    assert bu.check_version(bu.Version(1,0,0), bu.Version(2,0,0), bu.Version(bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT)) == False
    assert bu.check_version(bu.Version(1,0), bu.Version(2,0), bu.Version(bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT)) == False
    assert bu.check_version(bu.Version(1), bu.Version(2), bu.Version(bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT)) == False

    # unmeeting major conditions without min version
    assert bu.check_version(bu.Version(2,0,0), bu.Version(bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT), bu.Version(1,0,0)) == False
    assert bu.check_version(bu.Version(2,0), bu.Version(bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT), bu.Version(1,0,0)) == False
    assert bu.check_version(bu.Version(2), bu.Version(bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT), bu.Version(1,0,0)) == False

    # unmeeting major conditions on lower bound
    assert bu.check_version(bu.Version(1,0,0), bu.Version(2,0,0), bu.Version(3,0,0)) == False
    assert bu.check_version(bu.Version(1,0), bu.Version(2,0), bu.Version(3,0)) == False
    assert bu.check_version(bu.Version(1), bu.Version(2), bu.Version(3)) == False
    
    # unmeeting major conditions on upper bound
    assert bu.check_version(bu.Version(3,0,0), bu.Version(1,0,0), bu.Version(2,0,0)) == False
    assert bu.check_version(bu.Version(3,0), bu.Version(1,1), bu.Version(2,0)) == False
    assert bu.check_version(bu.Version(3), bu.Version(1), bu.Version(2)) == False


def test_minor():
    # meeting conditions without max version
    assert bu.check_version(bu.Version(1,10,0), bu.Version(1,1,0), bu.Version(bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT)) == True
    assert bu.check_version(bu.Version(1,10), bu.Version(1,1,0), bu.Version(bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT)) == True

    # meeting conditions without min version
    assert bu.check_version(bu.Version(1,1,0), bu.Version(bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT), bu.Version(1,10,0)) == True
    assert bu.check_version(bu.Version(1,1), bu.Version(bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT), bu.Version(1,10,0)) == True

    # meeting conditions on lower bound
    assert bu.check_version(bu.Version(1,0,0), bu.Version(1,0,0), bu.Version(1,10,0)) == True
    assert bu.check_version(bu.Version(1,0), bu.Version(1,0), bu.Version(1,10)) == True
    
    # meeting conditions on upper bound
    assert bu.check_version(bu.Version(1,10,0), bu.Version(1,0,0), bu.Version(1,10,0)) == True
    assert bu.check_version(bu.Version(1,10), bu.Version(1,0), bu.Version(1,10)) == True
    
    # meeting conditions in between
    assert bu.check_version(bu.Version(1,5,0), bu.Version(1,0,0), bu.Version(1,10,0)) == True
    

    # unmeeting major conditions without max version
    assert bu.check_version(bu.Version(1,0,0), bu.Version(1,10,0), bu.Version(bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT)) == False
    assert bu.check_version(bu.Version(1,0), bu.Version(1,10), bu.Version(bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT)) == False

    # unmeeting major conditions without min version
    assert bu.check_version(bu.Version(1,10,0), bu.Version(bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT), bu.Version(1,0,0)) == False
    assert bu.check_version(bu.Version(1,10), bu.Version(bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT), bu.Version(1,0)) == False

    # unmeeting major conditions on lower bound
    assert bu.check_version(bu.Version(1,0,0), bu.Version(1,10,0), bu.Version(1,0,0)) == False
    assert bu.check_version(bu.Version(1,0), bu.Version(1,10), bu.Version(1,0)) == False
    
    # unmeeting major conditions on upper bound
    assert bu.check_version(bu.Version(1,10,0), bu.Version(1,0,0), bu.Version(1,5,0)) == False
    assert bu.check_version(bu.Version(1,10), bu.Version(1,0), bu.Version(1,5)) == False


def test_patch():
    # meeting conditions without max version
    assert bu.check_version(bu.Version(1,0,10), bu.Version(1,0,1), bu.Version(bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT)) == True

    # meeting conditions without min version
    assert bu.check_version(bu.Version(1,0,10), bu.Version(bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT), bu.Version(1,0,10)) == True

    # meeting conditions on lower bound
    assert bu.check_version(bu.Version(1,0,1), bu.Version(1,0,1), bu.Version(1,0,10)) == True
    
    # meeting conditions on upper bound
    assert bu.check_version(bu.Version(1,0,10), bu.Version(1,0,0), bu.Version(1,0,10)) == True
    
    # meeting conditions in between
    assert bu.check_version(bu.Version(1,0,5), bu.Version(1,0,1), bu.Version(1,0,10)) == True
    

    # unmeeting major conditions without max version
    assert bu.check_version(bu.Version(1,0,0), bu.Version(1,0,10), bu.Version(bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT,bu.VERSION_MAX_DIGIT)) == False

    # unmeeting major conditions without min version
    assert bu.check_version(bu.Version(1,0,10), bu.Version(bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT,bu.VERSION_MIN_DIGIT), bu.Version(1,0,0)) == False

    # unmeeting major conditions on lower bound
    assert bu.check_version(bu.Version(1,0,0), bu.Version(1,0,10), bu.Version(1,0,10)) == False
    
    # unmeeting major conditions on upper bound
    assert bu.check_version(bu.Version(1,0,10), bu.Version(1,0,0), bu.Version(1,0,5)) == False
