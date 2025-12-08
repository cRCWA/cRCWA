#!/bin/bash

echo ""
printf "\033[1m   c R C W A\n"
echo ""
echo "   Automatic test suite"
echo "   by Davide Bucci 2013-2025"
printf "\033[0m\n"

if [ -z "$1" ]
then
    executable="../../bin/crcwa"
else
    executable=$1
fi


printf " NOTE: the following tests will run $executable\n\n"



test_fail=0

cd simplefor
    if ./test_simplefor.sh $executable
    then
        ok=1
    else
        test_fail=1
    fi
cd ..

cd fortest
    if ./test_fortest.sh $executable
    then
        ok=1
    else
        test_fail=1
    fi
cd ..

cd lot_of_layers
    if ./test_lot_of_layers.sh $executable
    then
        ok=1
    else
        test_fail=1
    fi
cd ..

cd iftest
    if ./test_iftest.sh $executable
    then
        ok=1
    else
        test_fail=1
    fi
cd ..

cd loadforif
    if ./test_loadforif.sh $executable
    then
        ok=1
    else
        test_fail=1
    fi
cd ..

cd calc
    if ./test_calc.sh $executable
    then
        ok=1
    else
        test_fail=1
    fi
cd ..

cd simplemode
    if ./test_simplemode.sh $executable
    then
        ok=1
    else
        test_fail=1
    fi
cd ..

cd bloch
    if ./test_bloch.sh $executable
    then
        ok=1
    else
        test_fail=1
    fi
cd ..

cd 1D_test_angle
    if ./test_1D_test_angle.sh $executable
    then
        ok=1
    else
        test_fail=1
    fi
cd ..

cd symtest
    if ./test_symmetry.sh $executable
    then
        ok=1
    else
        test_fail=1
    fi
cd ..

cd parallel
    if ./test_parallel.sh $executable
    then
        ok=1
    else
        test_fail=1
    fi
cd ..

cd matdevtests
    if ./test_matdev.sh $executable
    then
        ok=1
    else
        test_fail=1
     fi
cd ..

cd multiplepropagation
    if ./test_multiplepropagation.sh $executable
    then
        ok=1
    else
        test_fail=1
     fi
cd ..
cd multiplepropagation2
    if ./test_multiplepropagation2.sh $executable
    then
        ok=1
    else
        test_fail=1
     fi
cd ..
cd monitortest
    if ./test_monitor.sh $executable
    then
        ok=1
    else
        test_fail=1
     fi
cd ..
cd testorientation
    if ./testorientation.sh $executable
    then
        ok=1
    else
        test_fail=1
     fi
cd ..
cd testfile
    if ./testfile.sh $executable
    then
        ok=1
    else
        test_fail=1
     fi
cd ..
cd dowhiletest
    if ./dowhiletest.sh $executable
    then
        ok=1
    else
        test_fail=1
     fi
cd ..
cd calculatedz
    if ./test_calculatedz.sh $executable
    then
        ok=1
    else
        test_fail=1
     fi
cd ..

echo
if test $test_fail != 0
then
    printf "\033[1m\x1b[31m SOME TESTS FAILED \033[0m\n"
    exit 1
else
    printf "\x1b[32mAll is Ok! :-)\033[0m\n" 
fi
