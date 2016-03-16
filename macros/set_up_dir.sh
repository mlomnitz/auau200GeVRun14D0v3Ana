#!bin/bash
mkdir production log err csh list
find -maxdepth 1 -user $USER -exec chmod g+rw {} \;
ln -s $PWD/StRoot/macros/submitPicoMixEventsMaker.sh .
ln -s $PWD/StRoot/macros/submitPicoMixEventsMaker.xml .