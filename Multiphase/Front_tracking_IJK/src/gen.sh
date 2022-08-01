#!/bin/bash
export prepro="python $TRIO_U_ROOT/bin/KSH/preprocessor.py"
for i in `ls *.P`
do
  echo processing $i
  b=`basename $i .P`
  if test "$REGENERE_TOUT" == "1"; then
      j=
  else
      j=`find . -name $b -newer $i`
  fi

  if test -z "$j"
  then
    echo Generating $b
    $prepro $i $b
    [ $? -ne 0 ] && exit -1

    $TRIO_U_ROOT/exec/astyle/bin/astyle --style=gnu --indent=spaces=2 --indent-cases --align-reference=type $b 1>/dev/null

  else
    echo Found file "$j" newer than "$i", skipping
  fi
done
