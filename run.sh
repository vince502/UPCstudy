#!/bin/bash

root -l -b -q 'EfficiencyPurity.C(0)' >& log_p0.txt &
root -l -b -q 'EfficiencyPurity.C(1)' >& log_p1.txt &
root -l -b -q 'EfficiencyPurity.C(2)' >& log_p2.txt &
root -l -b -q 'EfficiencyPurity.C(3)' >& log_p3.txt &
root -l -b -q 'EfficiencyPurity.C(4)' >& log_p4.txt &
root -l -b -q 'EfficiencyPurity.C(5)' >& log_p5.txt &
