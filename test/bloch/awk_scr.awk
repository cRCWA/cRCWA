function abs(x){ return ((x < 0.0) ? -x : x) } BEGIN {comp=0; norm=0;} FNR==NR{a[FNR]=$4; norm+=$4*$4; next} a[FNR]!=$4{comp+=abs($4-a[FNR])} END{ norm=sqrt(norm/FNR); if (comp>2e-4*norm){ exit 1 } }
