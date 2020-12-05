import os
import pandas as pd 
import argparse


def convertToPdb(dirname,filename):

    '''
    Convert the saved configurations in dirname to a pdb file with name filename. 
    '''

    data=pd.read_csv(
        os.path.join(dirname,"particles.dat") , 
        delim_whitespace=True   
    )

    #data=data.query("mask==1");

    with open(filename,"w") as f:
        
        for index, row in data.iterrows():
                f.write("HETATM{:5d} {:4} {:>3} {:1}{:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}    {:>2}\n".format(index,"C","VAL","A",int(row["time"]),row["x"],row["y"],row["z"],1.00,row["particle"],"C") )


        for iP, particleData in  data.groupby("particle"):
            particleData=particleData.sort_values(by="time")
            particleData=particleData.reset_index()

            for i in range( len(particleData)-1):
                f.write( "CONECT{:>5d}{:>5d}\n".format( particleData.loc[i,"index"] , particleData.loc[i+1,"index"] ) )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("configuration", help="Directory containing the configurations")
    parser.add_argument("pdb", help="Output pdb file")

    args=parser.parse_args()

    convertToPdb(args.configuration,args.pdb)