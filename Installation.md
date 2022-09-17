# Installation of software

The installation of software is tested on Ubuntu 20.04 version but it supposed to work on any linux machine. I have used r.avaflow 2.4 version of software. The code for this version is availbe within this reposatory [here](https://github.com/iamtekson/GLOF-simulation/tree/main/r.avaflow).The software installation guide is as below,

## Grass gis and R installation

The installation document is available here: [`r.avaflow/grass7.install.sh`](https://github.com/iamtekson/GLOF-simulation/tree/main/r.avaflow/grass.install.sh). The file can be exicuted as below,

> For me grass 8.2 version was installed.

```sh
sh grass.install.sh
```

## r.avaflow Installation

To install the r.avaflow package, simply type following command in terminal,

> Don't forgot to change the `url` parameter `/path/to/r.avaflow/directory`

```sh
# To start the grass interface
grass

# To install the extension
g.extension extension=r.avaflow url=/path/to/r.avaflow/directory
```

If everything run smooth, then it is installed successfully. You are now ready to use the software. 
