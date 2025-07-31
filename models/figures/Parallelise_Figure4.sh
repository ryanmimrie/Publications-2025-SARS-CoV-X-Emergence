#!/bin/bash

SESSION="parallel_32"

for i in {1..32}; do
  sed -e "s/file_prefix <- 1/file_prefix <- $i/" \
      -e "s/parallelisation <- 1/parallelisation <- 32/" \
      Dust_Script_Figure4.R > Dust_Script_Figure4_${i}.R
done

# --- Start a new detached session ---
tmux new-session -d -s "$SESSION"

# --- Manual Grid Construction (4x8 = 32 panes) ---
tmux split-window -v -t "$SESSION":0.0
tmux split-window -v -t "$SESSION":0.0
tmux split-window -v -t "$SESSION":0.1
tmux split-window -v -t "$SESSION":0.0
tmux split-window -v -t "$SESSION":0.4
tmux split-window -v -t "$SESSION":0.4
tmux split-window -v -t "$SESSION":0.6

tmux split-window -h -t "$SESSION":0.0
tmux split-window -h -t "$SESSION":0.0
tmux split-window -h -t "$SESSION":0.2

tmux split-window -h -t "$SESSION":0.4
tmux split-window -h -t "$SESSION":0.4
tmux split-window -h -t "$SESSION":0.6

tmux split-window -h -t "$SESSION":0.8
tmux split-window -h -t "$SESSION":0.8
tmux split-window -h -t "$SESSION":0.10

tmux split-window -h -t "$SESSION":0.12
tmux split-window -h -t "$SESSION":0.12
tmux split-window -h -t "$SESSION":0.14

tmux split-window -h -t "$SESSION":0.16
tmux split-window -h -t "$SESSION":0.16
tmux split-window -h -t "$SESSION":0.18

tmux split-window -h -t "$SESSION":0.20
tmux split-window -h -t "$SESSION":0.20
tmux split-window -h -t "$SESSION":0.22

tmux split-window -h -t "$SESSION":0.24
tmux split-window -h -t "$SESSION":0.24
tmux split-window -h -t "$SESSION":0.26

tmux split-window -h -t "$SESSION":0.28
tmux split-window -h -t "$SESSION":0.28
tmux split-window -h -t "$SESSION":0.30

# --- Launch scripts in each pane ---
# Change the script name in '' to run a different command.
# This expects a "_#" suffix before the file name extension starting at 0.

NUM_PANES=$(tmux list-panes -t "$SESSION":0 | wc -l)

for i in $(seq 0 $((NUM_PANES - 1))); do
    FILE_NUM=$((i + 1))
    tmux send-keys -t "$SESSION":0."$i" "conda activate Odin" C-m
    tmux send-keys -t "$SESSION":0."$i" "Rscript Dust_Script_Figure4_${FILE_NUM}.R" C-m
done

# --- Attach to session (needs to be last) ---
tmux attach -t "$SESSION"
