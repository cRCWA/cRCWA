#!/bin/bash

printf "1 -thread timing:  "

time_1="$( TIMEFORMAT='%lR';time ( ../../../afmm test1.fmm >null.bak ) 2>&1 1>/dev/null )" 
echo "$time_1"

printf "2 -threads timing: "
time_2="$( TIMEFORMAT='%lR';time ( ../../../afmm test2.fmm >null.bak ) 2>&1 1>/dev/null )"
echo "$time_2"

printf "3 -threads timing: "
time_3="$( TIMEFORMAT='%lR';time ( ../../../afmm test3.fmm >null.bak ) 2>&1 1>/dev/null )" 
echo "$time_3"

printf "4 -threads timing: "
time_4="$( TIMEFORMAT='%lR';time ( ../../../afmm test4.fmm >null.bak ) 2>&1 1>/dev/null )"
echo "$time_4"
 
printf "5 -threads timing: "
time_5="$( TIMEFORMAT='%lR';time ( ../../../afmm test5.fmm >null.bak ) 2>&1 1>/dev/null )"
echo "$time_5"

printf "6 -threads timing: "
time_6="$( TIMEFORMAT='%lR';time ( ../../../afmm test6.fmm >null.bak ) 2>&1 1>/dev/null )" 
echo "$time_6"

printf "7 -threads timing: "
time_7="$( TIMEFORMAT='%lR';time ( ../../../afmm test7.fmm >null.bak ) 2>&1 1>/dev/null )" 
echo "$time_7"

printf "8 -threads timing: "
time_8="$( TIMEFORMAT='%lR';time ( ../../../afmm test8.fmm >null.bak ) 2>&1 1>/dev/null )" 
echo "$time_8"

printf "12-threads timing: "
time_8="$( TIMEFORMAT='%lR';time ( ../../../afmm test12.fmm >null.bak ) 2>&1 1>/dev/null )" 
echo "$time_8"

printf "16-threads timing: "
time_8="$( TIMEFORMAT='%lR';time ( ../../../afmm test16.fmm >null.bak ) 2>&1 1>/dev/null )" 
echo "$time_8"



rm null.bak
rm *.idx
rm test_results.txt
