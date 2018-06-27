#!/bin/bash
name="$(sqlite3 $2/IRC/inputs.db "select name from mopac where id=$1")"
mopacl $2"/IRC/minf_"$name".mop" 2>/dev/null
mopacl $2"/IRC/minr_"$name".mop" 2>/dev/null

