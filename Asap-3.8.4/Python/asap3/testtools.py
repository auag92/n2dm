#!/usr/bin/env python
"""
A collection of tools used to implement the ASAP test suite.

Name: testtools.py

Description: A collection of tools used to make it easier to write the
ASAP test suite.

Usage: Used mainly to test if tests have failed or succeeded:

         ReportTest(name, value, expected, deviation)

       checks that the obtained value is closer to the expected value
       than deviation.  Name is the name for the test.

Example usage: Run it on the command line to test the tools themselved

Expected result: A number of successful and failed tests, self explanatory.
"""


import sys

class _Reporter:
    """Check if a test has failed, reports the result and remembers failed tests.

    A single instance, ReportTest, is created in this module and
    should be used.  It is a callable object with one additional
    method.
    
    Call with four arguments:
        name: the name of the test.
        value: the result of the test.
        expected: the expected result of the test.
        deviation: the maximally allowed deviation.

    The object has the additional method Summary(exit=0), printing a
    summary.  If the optional argument exit is non-zero, a failed test
    causes the script to exit with exit the given exit code.

    """
    def __init__(self):
        self.Clear()
        self.exitaftererror = False
        
    def Clear(self):
        self.failed = []
        self.success = 0

    def __call__(self, name, value, expected, deviation, silent=False):
        "Report the result of a test."
        if abs(value - expected) <= deviation:
            # Success
            if not silent:
                print "%s: SUCCESS" % name
            self.success += 1
        else:
            # Failure
            error = ("%s: FAILED, got %.8g expected %.8g (max deviation: %.8g)"
                     % (name, value, expected, deviation))
            print error
            self.failed.append(error)
            if self.exitaftererror:
                raise AssertionError, error

    def BoolTest(self, name, value, silent=False):
        "Report the result of a boolean test"
        if value:
            if not silent:
                 print "%s: SUCCESS" % name
            self.success += 1
        else:
            error = ("%s: FAILED." % (name,))
            print error
            self.failed.append(error)
            if self.exitaftererror:
                raise AssertionError, "Exiting after error"

    def ExitAfterError(self, ex=True):
        self.exitaftererror = ex
        
    def GetNumberOfErrors(self):
        "Return the number of failed tests."
        return len(self.failed)

    def Summary(self, exit=1):
        if self.failed:
            report =  "%d tests succeeded, %d TESTS FAILED." % (self.success,
                                                                len(self.failed))
            print report
            print ""
            print "Failed tests:"
            for failure in self.failed:
                print "   ", failure
            if exit:
                # TestAll.py catches the exception below, so clean up before
                # TestAll calls the next test.
                self.Clear()  
                raise AssertionError, report
        else:
            print "ALL TESTS SUCCEEDED (%d tests were performed)" % self.success

# The main test object:
ReportTest = _Reporter()


if __name__ == "__main__":
    # Test the test code.
    print "Testing the test code.  5 tests should succeed, 5 should fail."
    ReportTest("Expected success 1", 42.1, 42.1, 0.0)
    ReportTest("Expected success 2", 42.101, 42.1, 0.001)
    ReportTest("Expected success 3", 42.0995, 42.1, 0.001)
    ReportTest("Expected success 4", -0.001, -0.01, 0.01)
    ReportTest("Expected success 5", 0.001, -0.005, 0.01)
    ReportTest("Expected failure 1", 42.100001, 42.1, 0.0)
    ReportTest("Expected failure 2", 42.11, 42.1, 0.001)
    ReportTest("Expected failure 3", 42.09, 42.1, 0.001)
    ReportTest("Expected failure 4", -0.001, -0.01, 0.0001)
    ReportTest("Expected failure 5", 0.001, -0.01, 0.0001)
    print ""
    ReportTest.Summary(exit=0)
    print ""
    print "IMPORTANT: This self test succeeded if EXACTLY FIVE TESTS FAILED."
    print ""
    if len(ReportTest.failed) == 5:
        print "This test succeeded"
    else:
        print "THIS TEST FAILED !!!"
        print "5 tests should have failed, %d tests did fail" % (len(ReportTest.failed),)
        sys.exit(1)
        
