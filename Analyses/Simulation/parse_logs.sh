#!/bin/bash

for f in *.log; do
    tail -200 "$f" | awk '!/^   */' > metrics_"$f"
done

