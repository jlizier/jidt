#!/bin/sh

if [ "$#" -lt "2" ]; then
  echo "Usage:\n\t./benchmark.sh tag iterator\n\nExample:\n\t./benchmark.sh \"kNN kernel\" \"seq 1 10\"\n";
  exit
fi

# Very convenient functions to plot the results of this script
# function plot { gnuplot -e "plot '$1' using 1:2; pause -1" ; }
# function fitline { gnuplot -e "set fit quiet; f(x)=a*x+b; fit f(x) '$1' u 1:2 via a, b; plot '$1' u 1:2, f(x) t sprintf('f(x) = %.2fx + %.2f', a, b); pause -1" ;}

tag=$1
vals=$(eval "$2")
echo $vals

# Warm-up GPU before taking measurements
warmup=5
for n in $warmup ; do
  ./perftest 10 >/dev/null
done

for n in $vals ; do
  dur=$(./perftest $n | grep -i "$tag" | cut -f2 -d: | awk '{print $1}') ;
  echo "$n\t$dur" >> gpu_perftimes.txt ;
done

