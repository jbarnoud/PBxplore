#!/bin/bash

# Install the last version of pandoc on travis-CI if it is
# not already available in the cache
if [[ ! -f $HOME/.cabal/bin/pandoc ]]
then
    cabal update
    cabal install pandoc==1.15.2
else
    echo "Get pandoc from the cache"
fi
