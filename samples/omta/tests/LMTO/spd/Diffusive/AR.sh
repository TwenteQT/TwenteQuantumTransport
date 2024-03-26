rm G.dat


for Len in 10 
do

cd L"$Len"

for sup in 4 5 6
do

cd "$sup"x"$sup"

rm G.dat


for subdir in cf-*
do

cd $subdir

awk <andout >>../G.dat '/Grand total/{i=$3}END{print i}'

cd ..

done

awk -v len="$Len" -v sup="$sup" <G.dat >>../../G.dat '{i+=$1}END{print len,sup*sup,1/(i/NR),i/NR}'

cd ..

done

cd ..

done

