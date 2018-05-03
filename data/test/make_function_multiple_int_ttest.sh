#!/bin/bash

testfile=multiple_int_ttest.tab
echo -e "peptide\tgo\tint1\tint2\tint3\tint4\tint5\tint6" > $testfile
echo -e "A\tGO:0008152\t10\t20\t15\t30\t21\t30" >> $testfile
echo -e "B\tGO:0022610\t20\t30\t20\t3500\t2000\t2000" >> $testfile