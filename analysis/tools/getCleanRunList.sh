#!/bin/bash

cat file.list.15GeV.clean | cut -d'/' -f 13 | sort | uniq > run.list