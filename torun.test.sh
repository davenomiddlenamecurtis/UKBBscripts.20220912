echo Haha I ran at startup
echo Haha I ran at startup >startup.test.out
dx upload startup.test.out
dx terminate $DX_JOB_ID

# test line here


