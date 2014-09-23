for i in {3000..20000..1000}
do
Qmsub -h 1 drawlinesmap.py $(expr $i - 1000) $i
done
Qmsub -h 1 drawlinesmap.py 100 2000
