**Run**
>python src/otero_precipitation.py
>python src/tests.py



**Git Setup**
On the remote server:
>mkdir aedes_aegypti.git
>cd aedes_aegypti.git/
>git --bare init --shared

On the local:
>git clone exequiel@zotac.hopto.org:/home/exequiel/repositories/aedes_aegypti.git
>touch readme.md
>git add readme.md
>git commit -m "initial commit"
>git push origin master
Total 18 (delta 0), reused 0 (delta 0)
To exequiel@zotac.hopto.org:/home/exequiel/repositories/aedes_aegypti.git
 * [new branch]      master -> master

*Not sure if I should on the local, create a repo and push it to the origin (instead of cloning, changing and pushing). like in https://feeding.cloud.geek.nz/posts/setting-up-centralied-git-repository/*


**To clone the project**
>git clone exequiel@zotac.hopto.org:/home/exequiel/repositories/aedes_aegypti.git
or
>git clone ssh://exequiel@zotac.hopto.org:8080/home/exequiel/repositories/aedes_aegypti.git

**To create a new branch**
>git checkout -b development
Commit if you have something
>git push -u origin development

**To delete a branch**
>git push origin --delete development
>git branch -d development

**To create a tag**
>git checkout master
>git tag -a v0.1 -m "Version: 0.1"
>git tag -l
>git push origin v0.1

**To merge**
First merge master -> development
>git checkout development
>git merge master development
>git push origin development

Now merge development -> master
>git checkout master
>git merge master development
>git push origin master

 *the branch order in the merge command is not important*

**Profiling**
>python -m cProfile -s cumtime src/tests.py  > p.txt

**Usage**
To save:
>python src/tests.py save

To compare against previous results:
>python src/tests.py compare data/test/previous_results/2018-04-17__09_40_45.csv data/test/previous_results/2018-04-17__09_39_47.csv
or
>ls data/test/previous_results/*.csv |  python src/tests.py compare
