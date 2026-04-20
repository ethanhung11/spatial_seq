from typing import Literal
from time import time, ctime
from datetime import timedelta

def stopwatch(taskname: str = "", start=time(), mode: Literal[-1, 0, 1] = 0):
    if mode == 0:
        print(
            (
                "\n---------------------------------------------\n"
                f"begin {taskname} at {ctime(time())}, "
                f"(total: {timedelta(seconds=time()-start)})\n"
            )
        )
    elif mode == 1:
        print(f"finished {taskname} in {timedelta(seconds=time()-start)}\n")
    elif mode == -1:
        pass
    else:
        raise ValueError("valid modes only include -1, 0, or 1")
    return time()