# bash mass rename files --> run with
# ---bash test.sh---

cd ./data/Background_test/MT30
# for i in *.jpg ; do mv "$i" "${i//[thing to be replaced] /[thing to replace it with]}" ; done
for i in *.jpg ; do mv "$i" "${i//_ }" ; done
