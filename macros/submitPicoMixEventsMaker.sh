#########################################################################
# File Name: submitPicoD0AnaMaker.sh
# Created Time: Fri 08 May 2015 03:15:35 AM PDT
#########################################################################
#!/bin/bash
#umask 664
prod=$(date +"%F_%H_%M")_$USER
path=$PWD
echo $prod
mkdir $PWD/production/$prod
chmod g+rw $PWD/production/$prod 
mkdir $PWD/log/$prod
chmod g+rw $PWD/log/$prod
mkdir $PWD/err/$prod
chmod g+rw $PWD/err/$prod

star-submit-template -template submitPicoMixEventsMaker.xml -entities listOfFiles=$1,productionId=$prod,prodPath=$PWD

chmod g+rw schedTemplateExp.xml

#Last bug to fix, alitle bit of ugly fix
ls -1Apl *.dataset | grep -v /\$ | awk -v user=$USER '$3==user{print}{}' > temp_$USER.list
i=1
for F in $(cat temp_$USER.list) ; do
    if [ $i -eq 9 ] 
    then
	#echo $F
	chmod g+rw $F
	i=1
    else
	let i+=1
    fi
done
rm temp_$USER.list

ls -1Apl *.session.xml | grep -v /\$ | awk -v user=$USER '$3==user{print}{}' > temp_$USER.list
i=1
for F in $(cat temp_$USER.list) ; do
    if [ $i -eq 9 ] 
    then
	#echo $F
	chmod g+rw $F
	i=1
    else
	let i+=1
    fi
done
rm temp_$USER.list

cd report
ls -1Apl | grep -v /\$ | awk -v user=$USER '$3==user{print}{}' > temp_$USER.list
i=1
for F in $(cat temp_$USER.list) ; do
    if [ $i -eq 9 ] 
    then
	#echo $F
	chmod g+rw $F
	i=1
    else
	let i+=1
    fi
done
rm temp_$USER.list
cd ../csh

ls -1Apl *.list | grep -v /\$ | awk -v user=$USER '$3==user{print}{}' > temp_$USER.list
i=1
for F in $(cat temp_$USER.list) ; do
    if [ $i -eq 9 ] 
    then
	#echo $F
	chmod g+rw $F
	i=1
    else
	let i+=1
    fi
done
rm temp_$USER.list

cd ../