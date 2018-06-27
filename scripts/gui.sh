module avail &>ma
anaconda2=`awk '{for(i=1;i<=NF;i++) if($i~"anaconda2") print $i}' ma`
if [ ! -z "$anaconda2" ]; then
   module load $anaconda2  2>/dev/null
   python=`which python`
   echo "#!"$python > ${TSSCDS}/bin/GUI.py
   cat ${TSSCDS}/bin/gui0.py >>${TSSCDS}/bin/GUI.py
else
   echo "#!/usr/bin/python" > ${TSSCDS}/bin/GUI.py
   cat ${TSSCDS}/bin/gui0.py >>${TSSCDS}/bin/GUI.py
fi

chmod +x ${TSSCDS}/bin/GUI.py
GUI.py

