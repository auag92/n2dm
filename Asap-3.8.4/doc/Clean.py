#!/usr/bin/env python
import os
from os.path import join, walk


def f(arg, dir, files):
    for f in files:
        if f.endswith('~') or f.startswith('.'):
            os.remove(join(dir, f))

walk('.', f, None)
