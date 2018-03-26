for i in `seq 36.29 0.001 36.31`; do
	FHACK=$i ./extract_img ../data/vv.dat {49,50}000
	mv y_norm.asc yna_$i.asc
done
