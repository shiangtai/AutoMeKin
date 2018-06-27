#!/bin/bash
name="$(sqlite3 $2/inputs.db "select name from mopac where id=$1")"
mopacl $2"/TSs/"$name"_thermo.mop" 2>/dev/null
mopacl $2"/IRC/"$name"_ircf.mop" 2>/dev/null
mopacl $2"/IRC/"$name"_ircr.mop" 2>/dev/null
