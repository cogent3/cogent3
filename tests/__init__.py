#!/usr/bin/env python
sub_modules = ["test_draw", "test_phylo"]

for sub_module in sub_modules:
    exec(f"from {__name__} import {sub_module}")
