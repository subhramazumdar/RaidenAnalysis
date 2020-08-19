for z in snapshots/*.json
do
    echo $z    
    cp $z .
    

done

for z in *.json
do
    
    python3 graph_shortest.py $z
    rm $z
    
done
