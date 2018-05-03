#!/bin/bash

testfile=function_multiple_int.tab
echo -e "peptide\tgo\tint1\tint2\tint3" > $testfile
echo -e "A\tGO:0008152\t10\t20\t70" >> $testfile
echo -e "B\tGO:0022610\t20\t30\t30" >> $testfile