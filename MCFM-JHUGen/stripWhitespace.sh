
#!/bin/sh

function stripws(){
  dos2unix $1
  sed -Ei 's/[ \t]+$//' "$1"
}


for aDir in QCDLoop TensorReduction src
do
   pushd $aDir
   for sDir in $(ls ./)
   do
      if [ -d $sDir ];then
         pushd $sDir
         for f in $(ls ./ | grep -v ".so" | grep -v ".o" | grep -v ".a" | grep -v ".sh")
         do
            if [ -f $f ];then
               stripws $f
            fi
         done
         popd
      fi
   done
   popd
done
