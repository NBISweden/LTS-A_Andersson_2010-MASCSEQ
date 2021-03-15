#!/usr/bin/env python

import pandas as pd


def read_samples(f):
    df = pd.read_csv(f, index_col=0)
    return df.to_dict(orient="index")