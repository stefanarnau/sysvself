#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:31:30 2024

@author: plkn
"""

# Imports
import pandas as pd
import os

# Paths
path_logfiles = "/mnt/data_fast/sysvself/logfiles/"

# Subject list
subject_list = ["VP02", "VP03", "VP04", "VP05"]

s = 0

# Load data
fn = os.path.join(path_logfiles, subject_list[s] + "_logEdgy.txt")
df = pd.read_csv(fn, sep=" ", header=2)