#!/bin/bash

for m in {101..200}
do
for z in {1..6}
do
	fn_out=$z"k""$m.out"
        jb_out=$z"k""$m.o"
        fn_err=$z"k""$m.e"
		subCommand="Rscript runeuler"$z"seed""$m.R"
        bsubCommand="bsub -W 240 -n 1 -J $z"iBGep"$m -o $jb_out -e $fn_err -M 4000 -R \"rusage[mem=4000]\" \"$subCommand >> $fn_out\""
		command="echo $subCommand > $fn_out; $bsubCommand"
		echo "$command"
		eval "$command"
done
done

