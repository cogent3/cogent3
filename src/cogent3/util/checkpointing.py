#!/usr/bin/env python
import os
import pickle
import time


class Checkpointer:
    def __init__(self, filename, interval=None, noisy=True) -> None:
        if interval is None:
            interval = 1800
        self.filename = filename
        self.interval = interval
        self.last_time = time.time()
        self.noisy = noisy

    def available(self):
        return self.filename is not None and os.path.exists(self.filename)

    def load(self):
        assert self.filename is not None, "check .available() first"
        with open(self.filename, "rb") as f:
            obj = pickle.load(f)
        self.last_time = time.time()
        return obj

    def record(self, obj, msg=None, always=False) -> None:
        if self.filename is None:
            return
        now = time.time()
        elapsed = now - self.last_time
        if always or elapsed > self.interval:
            if self.noisy and msg is not None:
                pass
            with open(self.filename, "wb") as f:
                pickle.dump(obj, f)
            self.last_time = now
