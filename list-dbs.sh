#!/bin/sh

printf "Human Databases\n---------------\n";
cut -f 1,2,3,5,6 /data/jyin/Databases/gstools-db/gstools-db-config_human.txt

printf "\nMouse Databases\n---------------\n";
cut -f 1,2,3,5,6 /data/jyin/Databases/gstools-db/gstools-db-config_mouse.txt
