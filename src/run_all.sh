#!/usr/bin/env bash
doit -n 50 -f dodo_round0.py &&
doit -n 50 -f dodo_analyze.py &&
doit -n 50 -f dodo_results.py &&
doit -n 10 -f dodo_promiscuous.py
