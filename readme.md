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

**Profiling**
>python -m cProfile -s cumtime src/tests.py  > p.txt
