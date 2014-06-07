#!/usr/bin/env python
import os
import pandas as pd

meta_fn = os.path.expanduser("~/vervet_project/metadata/163_population_ucla_id_taxon.csv")
meta_df = pd.read_csv(meta_fn, index_col=0)

def determine_pop(ind_series):
    if ind_series["species"]=="Chlorocebus sabaeus":
        return 'sab'
    elif ind_series["species"]=="Chlorocebus aethiops aethiops":
        return 'aet'
    elif ind_series["species"]=="Chlorocebus tantalus":
        return 'tan'
    elif ind_series["species"]=="Chlorocebus pygerythrus cynosurus":
        return 'cyn'
    elif ind_series["species"]=="Chlorocebus pygerythrus pygerythrus":
        if ind_series["country"] in ["Kenya","Tanzania"]:
            return 'pyn'
        elif ind_series["country"] in ["South Africa","Botswana"]:
            return "pys"
        else:
            print "unknown pyg from", ind_series["country"]
    else:
        print "Could not infer population for", ind_series

cols = meta_df.columns.values
cols[0] = "country"
cols[1] = "species"
meta_df.columns = cols        
        
meta_df["population"] = meta_df.apply(determine_pop,axis=1)

meta_df.to_csv(meta_fn)

id_fn_base = os.path.expanduser("~/vervet_project/metadata/163_pop_ucla_id_")

for pop,df in meta_df.groupby("population"):
    fn = id_fn_base + pop + ".txt"
    with open(fn,"w") as f:
        for id in df.index.values:
            f.write(id + "\n")

