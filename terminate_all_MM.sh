#!/bin/bash
#
# Terminate grid jobs
#
#
# mikael.mieskolainen@cern.ch, 2018


rm *.jdl *.root *.xml

## DATA
#aliroot -b -q .x runGridMM.C'("terminate", "LHC15f", "pass2")'
#aliroot -b -q .x runGridMM.C'("terminate", "LHC15h", "pass1")'
#aliroot -b -q .x runGridMM.C'("terminate", "LHC17j", "pass1")'


## MC 2015
#aliroot -b -q .x runGridMM.C'("terminate", "LHC15g3a2")'
#aliroot -b -q .x runGridMM.C'("terminate", "LHC15g3a3")'
#aliroot -b -q .x runGridMM.C'("terminate", "LHC15g3c2")'


## MC 2016
#aliroot -b -q .x runGridMM.C'("terminate", "LHC16a2a2")'
#aliroot -b -q .x runGridMM.C'("terminate", "LHC16a2b2")'
#aliroot -b -q .x runGridMM.C'("terminate", "LHC16a2c2")'
#aliroot -b -q .x runGridMM.C'("terminate", "LHC16a2d2")'
#aliroot -b -q .x runGridMM.C'("terminate", "LHC16a2d2_plus")'


## MC 2017
aliroot -b -q .x runGridMM.C'("terminate", "LHC17h7a")'
aliroot -b -q .x runGridMM.C'("terminate", "LHC17h7b")'

