## Data Driven Reachability Analysis using Matrix Zonotopes
<br /> 
This repo cotains the code for our paper:<br /> 
[1] Amr Alanwar, Anne Koch, Frank Allgöwer, Karl Johansson "Data Driven Reachability Analysis Using Matrix Zonotopes <br />



We consider the problem of reachability analysis from noisy data given that the system 
model is unknown. Identifying a model is a preliminary step in the state of the art 
reachability analysis approaches. However, systems are becoming more complex and data 
is becoming more readily available.<br />

<br /> <br />
<p align="center">
<img
src="Figures/idea.png"
raw=true
alt="Subject Pronouns"
width=500
/>
</p>

## Running 
1- Download [CORA](https://github.com/TUMcps/CORA) and [MPT](https://www.mpt3.org) toolboxs.<br />
2- Add CORA nad MPT folder and subfolders to the Matlab path.  <br />
3- Add the repo folder and subfolders to the Matlab path.  <br />
4- run t_linearDT.m for linear system.<br />
5- run t_nonlinearDT.m for nonlinear system.<br />
<br />

Our paper Bibtex is as follows:<br />
```
@article{datadriven_reach,
  title={Data-Driven Reachability Analysis Using Matrix Zonotopes},
  author={Alanwar, Amr and Koch, Anne and Allgöwer, Frank and Johansson, Karl Henrik},
  journal={arXiv preprint arXiv:2011.08472},
  year={2020}
}
```
