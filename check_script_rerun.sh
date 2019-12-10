#! /bin/bash
#CMD="source /Users/majoburo/.bash_profile; tmux new-session -d -s Diagnosis \"source /Users/majoburo/.bash_profile; cd work/Diagnosis/; python diagnose.py recipients.py \" "
CMD="source /Users/majoburo/.zshrc; cd /Users/majoburo/work/Diagnosis/; nohup python diagnose.py recipients.py 2>&1 >> diagnose.stdout & "
SEARCHSTR="python diagnose.py" # should be start of last part of command
echo $CMD

# run pgrep to see if command exists
pidcmd="pgrep -f \"${SEARCHSTR}\" "
eval $pidcmd
status=$?
if [ $status -ne 0 ] ; then
    echo "Script dead, rerunning"
    echo $CMD
    eval $CMD
else
    echo "Running."
fi
