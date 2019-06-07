#!/usr/bin/env python

import os, sys, time
sys.path.append(os.pardir)

import check_obsysrc
import filecmp
import shutil
import unittest

infile = os.path.join("input","obsys-test-gaas.rc")
check_obsysrc = os.path.join(os.pardir, "check_obsysrc.py")

class CheckObsysRcTest(unittest.TestCase):

    #.......................................................................
    def test_1(self):
        "Test check_obsysrc without any flags"
        global check_obsysrc, infile

        newfil = os.path.join("outdir", "obsys-test1.rc.new")
        newexp = os.path.join("outexp", "obsys-test1.rc.new")

        errfil  = os.path.join("outdir", "obsys-test1.rc.err")
        errfil_ = os.path.join("outdir", "obsys-test1.rc.err_")
        errexp_ = os.path.join("outexp", "obsys-test1.rc.err_")

        if os.path.isfile(newfil): os.remove(newfil)
        if os.path.isfile(errfil): os.remove(errfil)

        cmd = check_obsysrc+" "+infile+" "+newfil+" "+errfil
        print(cmd)
        os.system(cmd)

        with open(errfil, mode="r") as input:
            with open(errfil_, mode="w") as output:
                for line in input:
                    if line.find("gap start") == -1:
                        output.write(line)

        self.assertTrue(filecmp.cmp(newfil, newexp))
        self.assertTrue(filecmp.cmp(errfil_, errexp_))

        os.remove(newfil)
        os.remove(errfil)
        os.remove(errfil_)

    #.......................................................................
    def test_2(self):
        'Test with "--ignore_gaps mod04_006_flk" flag'
        global check_obsysrc, infile

        newfil = os.path.join("outdir", "obsys-test2.rc.new")
        newexp = os.path.join("outexp", "obsys-test2.rc.new")

        errfil  = os.path.join("outdir", "obsys-test2.rc.err")
        errfil_ = os.path.join("outdir", "obsys-test2.rc.err_")
        errexp_ = os.path.join("outexp", "obsys-test2.rc.err_")

        if os.path.isfile(newfil): os.remove(newfil)
        if os.path.isfile(errfil): os.remove(errfil)

        flags = " --ignore_gaps mod04_006_flk "
        cmd = check_obsysrc+flags+infile+" "+newfil+" "+errfil
        print(cmd)
        os.system(cmd)

        with open(errfil, mode="r") as input:
            with open(errfil_, mode="w") as output:
                for line in input:
                    if line.find("gap start") == -1:
                        output.write(line)

        self.assertTrue(filecmp.cmp(newfil, newexp))
        self.assertTrue(filecmp.cmp(errfil_, errexp_))

        os.remove(newfil)
        os.remove(errfil)
        os.remove(errfil_)

    #.......................................................................
    def test_3(self):
        'Test with "--obslist mod04_006_flk --ignore_gaps all=2" flags'
        global check_obsysrc, infile

        newfil = os.path.join("outdir", "obsys-test3.rc.new")
        newexp = os.path.join("outexp", "obsys-test3.rc.new")

        errfil  = os.path.join("outdir", "obsys-test3.rc.err")
        errfil_ = os.path.join("outdir", "obsys-test3.rc.err_")
        errexp_ = os.path.join("outexp", "obsys-test3.rc.err_")

        if os.path.isfile(newfil): os.remove(newfil)
        if os.path.isfile(errfil): os.remove(errfil)

        flags = " --obslist mod04_006_flk --ignore_gaps all=2 "
        cmd = check_obsysrc+flags+infile+" "+newfil+" "+errfil
        print(cmd)
        os.system(cmd)

        with open(errfil, mode="r") as input:
            with open(errfil_, mode="w") as output:
                for line in input:
                    if line.find("gap start") == -1:
                        output.write(line)

        self.assertTrue(filecmp.cmp(newfil, newexp))
        self.assertTrue(filecmp.cmp(errfil_, errexp_))

        os.remove(newfil)
        os.remove(errfil)
        os.remove(errfil_)

    #.......................................................................
    def test_4(self):
        'Test with "--obslist mod04_006_flk --ignore_gaps all=0" flags'
        global check_obsysrc, infile

        newfil = os.path.join("outdir", "obsys-test4.rc.new")
        newexp = os.path.join("outexp", "obsys-test4.rc.new")

        errfil  = os.path.join("outdir", "obsys-test4.rc.err")
        errfil_ = os.path.join("outdir", "obsys-test4.rc.err_")
        errexp_ = os.path.join("outexp", "obsys-test4.rc.err_")

        if os.path.isfile(newfil): os.remove(newfil)
        if os.path.isfile(errfil): os.remove(errfil)

        flags = " --obslist mod04_006_flk --ignore_gaps all=0 "
        cmd = check_obsysrc+flags+infile+" "+newfil+" "+errfil
        print(cmd)
        os.system(cmd)

        with open(errfil, mode="r") as input:
            with open(errfil_, mode="w") as output:
                for line in input:
                    if line.find("gap start") == -1:
                        output.write(line)

        self.assertTrue(filecmp.cmp(newfil, newexp))
        self.assertTrue(filecmp.cmp(errfil_, errexp_))

        os.remove(newfil)
        os.remove(errfil)
        os.remove(errfil_)

    #.......................................................................
    def test_5(self):
        "Test should fail if newfil equals input file"
        global check_obsysrc, infile

        errfil  = os.path.join("outdir", "obsys-test5.rc.err")
        self.assertRaises(Exception, check_obsysrc, infile, infile, errfil)

    #.......................................................................
    def test_6(self):
        'Test with "--obslist nim07_tomseff_nc" flag'
        global check_obsysrc, infile

        newfil = os.path.join("outdir", "obsys-test6.rc.new")
        newexp = os.path.join("outexp", "obsys-test6.rc.new")

        errfil  = os.path.join("outdir", "obsys-test6.rc.err")
        errfil_ = os.path.join("outdir", "obsys-test6.rc.err_")
        errexp_ = os.path.join("outexp", "obsys-test6.rc.err_")

        if os.path.isfile(newfil): os.remove(newfil)
        if os.path.isfile(errfil): os.remove(errfil)

        flags = " --obslist nim07_tomseff_nc " 
        cmd = check_obsysrc+flags+infile+" "+newfil+" "+errfil
        print(cmd)
        os.system(cmd)

        with open(errfil, mode="r") as input:
            with open(errfil_, mode="w") as output:
                for line in input:
                    if line.find("gap start") == -1:
                        output.write(line)

        self.assertTrue(filecmp.cmp(newfil, newexp))
        self.assertTrue(filecmp.cmp(errfil_, errexp_))

        os.remove(newfil)
        os.remove(errfil)
        os.remove(errfil_)

#.......................................................................
if __name__ == "__main__":

    # setup once for all tests
    #-------------------------
    print("running setup")

    if os.path.isdir("outdir"):
        print("> removing outdir directory")
        shutil.rmtree("outdir")
    print("> making outdir directory")
    os.mkdir("outdir")

    if os.path.isdir("archive"):
        print("> removing archive directory")
        shutil.rmtree("archive")

    tarfile = os.path.join("input", "archive.tar")
    print("> untarring "+tarfile)
    os.system("tar xf "+tarfile)

    # run tests
    #----------
    unittest.main(exit=False)

    # clean up afterwards
    #--------------------
    if os.path.isdir("archive"):
        print("removing archive directory")
        shutil.rmtree("archive")
