#!/bin/bash
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

ssh -fN -L 2222:pchuckle:22 dcurtis@tails.cs.ucl.ac.uk
sftp -P2222 rejudcu@localhost

# must kill ssh later



