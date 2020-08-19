for z in snapshots/*.json
do
    echo $z    
    cp $z .
    

done

for z in *.json
do
    
    python3 degree_dist.py $z output.csv
    rm $z
    
done
