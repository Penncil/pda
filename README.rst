================
pda installation 
================

Build package::  
 
 cd build
 ./mkpda.sh

=======================
pda docker installation 
=======================

Step 1: install docker and docker-compose https://docs.docker.com/compose/install/

Step 2: Build container::  
 
 docker-compose build
    
Step 3: Launch containers::

 docker-compose up

Step 4: Launch R::

 docker-compose exec site R
