#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd

data = pd.read_csv("25GPM-25PSIG.csv", sep="\t")
time = data.loc[:,"Time (sec)"]

