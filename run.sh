#!/bin/bash
echo 'compiling...'
set -xe
libs=`echo .:./lib/lwjgl/*:./lib/lwjgl-opengl/*:./lib/lwjgl-glfw/*:./lib/joml-1.9.6.jar`
classfiles=`echo .:./shader/*:./world/*:./window/*`
echo ${libs}:${classfiles} > .classpath

if [[ "$1" =~ .java ]]; then
    javac -cp ${libs}:${classfiles} $1
    exit $?
fi
javac -cp ${libs}:${classfiles} */*.java
javac -cp ${libs}:${classfiles} Main.java
if [ 0 != $? ]; then
    exit 1
fi
echo 'running...'
java -cp .:${libs}:${classfiles} Main
