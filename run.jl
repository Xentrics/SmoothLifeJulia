#!/usr/bin/env julia

## Tested with Julia

include("main_functions.jl")

res = simulate(runs=100, r_in = 8, r_out = 20, dt = 0.05, sleepTime = 0);

save_runs_as_video(res["all"], "all_steps_test", ".mp4")
