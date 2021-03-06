#!/usr/bin/python
import sys
import os
import commands
from commands import getstatusoutput
import datetime
import argparse
import string
if __name__ == '__main__':
    parser = argparse.ArgumentParser (description = 'produce ntuples with WW semileptonic final state')
    parser.add_argument ('-i', '--inputFolder' , default = '../AllHadronicSUSY/TreeMaker/' , help='input folder with the reduced trees')
    parser.add_argument ('-o', '--output' , default = 'RSGraviton4000.root', help='output file')
    parser.add_argument ('-mc', '--ismc' , default = '0', help='is MC or not')
    parser.add_argument ('-l', '--lepton' , default = 'mu', help='lepton category (mu or el)')
    parser.add_argument ('-t', '--tree' , default = 'TreeMaker2/PreSelection', help='name of the input tree')
    parser.add_argument ('-n', '--name' , default = 'ReducedSelection_RSGraviton4000.root' , help='input file')
    parser.add_argument ('-w', '--xsecWeight' , default = '0.0002739' , help='xsec (pb)')
    parser.add_argument ('-no', '--numberOfEntries' , default = '25000' , help='number of initial entries of the dataset')
    parser.add_argument ('-a', '--aqgc', default = '0', help='includes AQGC reweigh factors')
    parser.add_argument ('-f', '--firstEvent' , default = '1' , help='first event to process')
    parser.add_argument ('-j', '--PUJetID' , default = '0' , help='Apply PUJetID to AK4CHS jets')
    args = parser.parse_args ()

    print 'produceWWNtuples '+args.inputFolder+' '+args.output+' '+args.ismc+' '+args.lepton+' '+args.tree+' '+args.name+' '+args.xsecWeight+' '+args.numberOfEntries+' '+args.aqgc+' '+args.firstEvent+' '+args.PUJetID
    os.system('./produceWWNtuples.exe '+args.inputFolder+' '+args.output+' '+args.ismc+' '+args.lepton+' '+args.tree+' '+args.name+' '+args.xsecWeight+' '+args.numberOfEntries+' '+args.aqgc+' '+args.firstEvent+' '+args.PUJetID)
