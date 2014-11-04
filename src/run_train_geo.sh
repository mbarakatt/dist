for i in {0..2500..50} 
do
	Qmsub -h 3 -n 1 train_geo.py $i $(expr $i + 50 )
done

