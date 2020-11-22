
txval=10000
budget=3000000

while penalty<=0.1:
    for z in snapshots/*.json:
            python3 graph_main_attack.py z budget penalty txval output_synthetic_test.csv
            
      

