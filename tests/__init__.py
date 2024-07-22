#!/usr/bin/env python
import os
import pathlib

os.chdir(pathlib.Path(__file__).parent)

sub_modules = ["test_draw", "test_phylo"]

for sub_module in sub_modules:
    exec(f"from {__name__} import {sub_module}")
