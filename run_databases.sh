#!/bin/bash

sbatch /y/people/jsd523/dev/privateer/update_privateer_database.job
date=$(date +"%Y-%m-%d %T")

curl -X POST "https://chat.googleapis.com/v1/spaces/AAAANjHPkNs/messages?key=AIzaSyDdI0hCZtE6vySjMm-WEfRq3CPzqKqqsHI&token=cPXLED9TFEasmCAmoSab9D0F1AqPc5rV37haAb1GLEY" -H "Content-Type: application/json" -d "{'text': 'Queued Database Jobs - `date`'}"