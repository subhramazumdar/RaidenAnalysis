for z in snapshots/*.json
do
    echo $z    
    cp $z .
    

done

for z in *.json
do
    
    python3 graph_capacityloss.py $z
    rm $z
    
done
