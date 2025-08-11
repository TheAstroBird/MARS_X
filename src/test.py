import pandas as pd
import matplotlib.pyplot as plt
import subprocess as sub

for composition in ['BF97', 'MA79', 'LF97', 'S99', 'KC08', 'T13']:
    for areoterm in ['ATL', 'ATM', 'ATH']:
        for i in range(3):
            inp = composition + '\n' + areoterm + '\n' + str(i+1) + '\n'
            sub.run('python3 perple2xlsx.py', input=inp, shell=True, text=True)