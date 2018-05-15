#!/bin/bash

testfile=function_cog_multiple_int_ttest_filtering.tab
echo -e "peptide\tcog\tint1\tint2\tint3\tint4\tint5\tint6" > $testfile
echo -e "A\tO\t0\t0\t15\t30\t21\t30" >> $testfile
echo -e "B\tC\t20\t30\t20\t3500\t2000\t2000" >> $testfile